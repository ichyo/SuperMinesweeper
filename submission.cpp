#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <string>
#include <cassert>
#include <cstring>
#include <tuple>
#include <algorithm>
#include <chrono>
#include <bitset>
#include <functional>

#ifdef LOCAL_TEST
#include "./dbg.hpp"
#define MAX_RUNTIME 7500
#else
#define dbg(...)
#define MAX_RUNTIME 10000
#endif

#define REP(i,n) for(int i = 0; i < (n); ++i)

using namespace std;

const int MAX_N = 50;

const int MAX_POINTS = 64 * 6;

const int ESTIMATED_ADD_BOMB = 4; // TODO: tuning here

using Pos = pair<int, int>;

struct Constraint {
    int value;
    set<Pos> positions;

    Constraint(int value_, set<Pos> positions_): value(value_), positions(positions_) {}

    int max_one() const {
        return value;
    }
    int max_zero() const {
        return positions.size() - value;
    }
};

using Constraints = vector<Constraint>;

// https://github.com/atcoder/ac-library/blob/master/atcoder/dsu.hpp
// CC0 license
struct UnionFind {
    public:
    UnionFind() : _n(0) {}
    UnionFind(int n) : _n(n), parent_or_size(n, -1) {}

    int merge(int a, int b) {
        assert(0 <= a && a < _n);
        assert(0 <= b && b < _n);
        int x = leader(a), y = leader(b);
        if (x == y)
            return x;
        if (-parent_or_size[x] < -parent_or_size[y])
            std::swap(x, y);
        parent_or_size[x] += parent_or_size[y];
        parent_or_size[y] = x;
        return x;
    }

    bool same(int a, int b) {
        assert(0 <= a && a < _n);
        assert(0 <= b && b < _n);
        return leader(a) == leader(b);
    }

    int leader(int a) {
        assert(0 <= a && a < _n);
        if (parent_or_size[a] < 0)
            return a;
        return parent_or_size[a] = leader(parent_or_size[a]);
    }

    int size(int a) {
        assert(0 <= a && a < _n);
        return -parent_or_size[leader(a)];
    }

    std::vector<std::vector<int>> groups() {
        std::vector<int> leader_buf(_n), group_size(_n);
        for (int i = 0; i < _n; i++) {
            leader_buf[i] = leader(i);
            group_size[leader_buf[i]]++;
        }
        std::vector<std::vector<int>> result(_n);
        for (int i = 0; i < _n; i++) {
            result[i].reserve(group_size[i]);
        }
        for (int i = 0; i < _n; i++) {
            result[leader_buf[i]].push_back(i);
        }
        result.erase(
                std::remove_if(result.begin(), result.end(),
                    [&](const std::vector<int> &v) { return v.empty(); }),
                result.end());
        return result;
    }

    private:
    int _n;
    // root node: -1 * component size
    // otherwise: parent
    std::vector<int> parent_or_size;
};

struct Params {
    int N;
    int M;
    int D;
};

class Cell {
    const int HIDDEN = -1;
    const int BOMB = -2;
    const int SAFE = -3;

    int value;

    public:
    Cell() {
        value = HIDDEN;
    }

    bool is_hidden() const {
        return value == HIDDEN;
    }

    bool is_bomb() const {
        return value == BOMB;
    }

    bool has_value(int x) {
        assert(x >= 0);
        return value == x;
    }

    bool is_value() const {
        return value >= 0;
    }

    bool is_safe() const {
        return value == SAFE;
    }

    int get_value() const {
        assert(is_value());
        return value;
    }

    void set_value(int value) {
        assert(is_hidden() || is_safe());
        assert(value >= 0);
        this->value = value;
    }

    void set_bomb() {
        assert(is_hidden());
        this->value = BOMB;
    }

    void set_safe() {
        assert(is_hidden());
        this->value = SAFE;
    }
};

class Command {
    enum Type {
        Stop,
        Open,
        Flag,
    };

    Type type;
    int _row;
    int _col;

    Command(Type type) {
        this->type = type;
    }

    Command(Type type, int row, int col) {
        this->type = type;
        this->_row = row;
        this->_col = col;
    }

    public:
    static Command stop() {
        return Command(Type::Stop);
    }

    static Command open(int r, int c) {
        return Command(Type::Open, r, c);
    }

    static Command flag(int r, int c) {
        return Command(Type::Flag, r, c);
    }

    int row() const {
        return _row;
    }

    int col() const {
        return _col;
    }

    bool is_stop() {
        return type == Type::Stop;
    }

    bool is_open() {
        return type == Type::Open;
    }

    bool is_flag() {
        return type == Type::Flag;
    }

    string to_string() {
        if (is_stop()) {
            return "STOP";
        }
        if (is_open()) {
            return "G " + ::to_string(_row) + " " + ::to_string(_col);
        }
        if (is_flag()) {
            return "F " + ::to_string(_row) + " " + ::to_string(_col);
        }
        assert(false);
    }
};

struct Stats {
    int guess_count = 0;
    int random_guess_count = 0;

    void print() {
        cerr << "guess_count: " << guess_count << endl;
        cerr << "random_guess_count: " << random_guess_count << endl;
    }
};

class Solver {
    Params params;

    Stats stats;

    int grid_unknown_count = 0;
    int grid_bomb_count = 0;
    int score_uncover = 0;
    int score_mine_hit = 0;

    bool backtrace_timeout = false;

    vector<vector<Cell>> judge_grid;
    vector<vector<Cell>> grid;

    long runtime = 0;
    std::chrono::system_clock::time_point last_update = chrono::system_clock::now();
    std::chrono::system_clock::time_point now = chrono::system_clock::now();
    int now_counter = 0;

    void set_runtime(long runtime) {
        this->runtime = runtime;
        last_update = chrono::system_clock::now();
        now = last_update;
        now_counter = 0;
    }

    long get_runtime(bool force=false) {
        if (force || ++now_counter == 1000) {
            now_counter = 0;
            now = chrono::system_clock::now();
        }
        return runtime + chrono::duration_cast<std::chrono::milliseconds>(now-last_update).count();
    }

    vector<Command> next_commands;

    int count_safe_neighbors[MAX_N][MAX_N];
    int count_bomb_neighbors[MAX_N][MAX_N];
    set<Pos> unknown_neighbors[MAX_N][MAX_N];
    vector<Pos> neighbors[MAX_N][MAX_N];

    set<Pos> positions_with_unknown_values;

    bool with_in(int r1, int c1, int r2, int c2, int D) {
        return (r1 - r2) * (r1 - r2) + (c1 - c2) * (c1 - c2) <= D;
    }

    double calc_uncover_ratio() {
        const int all = params.N * params.N - params.M;
        return (double)score_uncover/(double)all;
    }

    double score() {
        return calc_uncover_ratio() / (1.0 + score_mine_hit);
    }

    public:
    Solver(int N, int M, int D, int row, int col) {
        params.N = N;
        params.M = M;
        params.D = D;

        grid = vector<vector<Cell>>(N, vector<Cell>(N, Cell()));
        judge_grid = vector<vector<Cell>>(N, vector<Cell>(N, Cell()));

        memset(count_safe_neighbors, 0, sizeof(count_safe_neighbors));
        memset(count_bomb_neighbors, 0, sizeof(count_bomb_neighbors));
        REP(r, N) REP(c, N) {
            REP(nr, N) REP(nc, N) {
                if (with_in(r, c, nr, nc, D)) {
                    neighbors[r][c].push_back(make_pair(nr, nc));
                    unknown_neighbors[r][c].insert(make_pair(nr, nc));
                }
            }
        }

        grid_unknown_count = params.N * params.N;

        judge_update_value(row, col, 0, 0);
    }

    bool is_game_end() {
        const int all = params.N * params.N - params.M;
        return all == score_uncover;
    }

    bool is_last_move() {
        const int all = params.N * params.N - params.M;
        return all == score_uncover + 1;
    }

    Command choose_random_unknown() {
        while(true) {
            int r = rand() % params.N;
            int c = rand() % params.N;
            if (grid[r][c].is_hidden() && judge_grid[r][c].is_hidden()) {
                return Command::open(r, c);
            }
        }
    }

    Constraints extract_constraints() {
        Constraints results;
        for(auto p : positions_with_unknown_values) {
            int r, c;
            tie(r, c) = p;
            assert(grid[r][c].is_value() && !unknown_neighbors[r][c].empty());
            const int value = grid[r][c].get_value() - count_bomb_neighbors[r][c];
            assert(unknown_neighbors[r][c].size() > value);
            results.push_back(Constraint(value, set<Pos>(unknown_neighbors[r][c].begin(), unknown_neighbors[r][c].end())));
        }
        return results;
    }

    vector<Constraints> split_constraints(const Constraints& constraints) {
        UnionFind uf(constraints.size());
        vector<vector<int>> c_id(params.N, vector<int>(params.N, -1));
        for(int i = 0; i < constraints.size(); i++) {
            for (const auto& p : constraints[i].positions) {
                const int j = c_id[p.first][p.second];
                if (j == -1) {
                    c_id[p.first][p.second] = i;
                } else {
                    uf.merge(i, j);
                }
            }
        }
        const auto groups = uf.groups();
        vector<Constraints> result;
        for(const auto& group : groups) {
            Constraints list;
            for(int i : group) {
                list.push_back(constraints[i]);
            }
            result.push_back(list);
        }
        return result;
    }

    void dfs(
        int index,
        bitset<MAX_POINTS>& is_mine,
        vector<int>& count_zero,
        vector<int>& count_one,
        vector<bitset<MAX_POINTS>>& results,
        const Constraints& constraints,
        const vector<vector<int>>& p2c
    ) {
        if (index == p2c.size()) {
            results.push_back(is_mine);
            return;
        }

        if (backtrace_timeout) {
            return;
        }

        if (get_runtime() > MAX_RUNTIME * 0.95) {
            dbg(get_runtime());
            dbg(p2c.size());
            dbg(constraints.size());
            dbg(results.size());
            backtrace_timeout = true;
            return;
        }

        {
            int violate = p2c[index].size();
            for (int i = 0; i < p2c[index].size(); i++) {
                const int x = p2c[index][i];
                if (++count_zero[x] > constraints[x].max_zero()) {
                    violate = i;
                    break;
                }
            }
            if (violate == p2c[index].size()) {
                is_mine[index] = false;
                dfs(index + 1, is_mine, count_zero, count_one, results, constraints, p2c);
                violate--;
            }
            for (int i = 0; i <= violate; i++) {
                const int x = p2c[index][i];
                --count_zero[x];
            }
        }
        {
            int violate = p2c[index].size();
            for (int i = 0; i < p2c[index].size(); i++) {
                const int x = p2c[index][i];
                if (++count_one[x] > constraints[x].max_one()) {
                    violate = i;
                    break;
                }
            }
            if (violate == p2c[index].size()) {
                is_mine[index] = true;
                dfs(index + 1, is_mine, count_zero, count_one, results, constraints, p2c);
                violate--;
            }
            for (int i = 0; i <= violate; i++) {
                const int x = p2c[index][i];
                --count_one[x];
            }
        }
    }

    void backtrace(const Constraints& constraints, vector<pair<double, Pos>>& probabilities) {
        if (backtrace_timeout) {
            return;
        }

        set<Pos> point_set;
        for(const auto& c : constraints) {
            for(const auto& p : c.positions) {
                point_set.insert(p);
            }
        }
        vector<Pos> points(point_set.begin(), point_set.end());
        if (points.size() > MAX_POINTS) {
            return;
        }

        vector<vector<int>> p2c(points.size(), vector<int>());
        for(int i = 0; i < constraints.size(); i++) {
            for (int j = 0; j < points.size(); j++) {
                if (constraints[i].positions.count(points[j])) {
                    p2c[j].push_back(i);
                }
            }
        }
        bitset<MAX_POINTS> is_mine;
        vector<int> count_zero(constraints.size(), 0);
        vector<int> count_one(constraints.size(), 0);

        vector<bitset<MAX_POINTS>> results;
        dfs(0, is_mine, count_zero, count_one, results, constraints, p2c);

        if (backtrace_timeout) {
            return;
        }

        assert(!results.empty());

        for (int i = 0; i < points.size(); i++) {
            const auto& p = points[i];
            int r, c;
            tie(r, c) = p;

            if (!grid[r][c].is_hidden()) {
                continue;
            }

            int count_one = 0;
            for(const auto& result : results) {
                if (get_runtime() > MAX_RUNTIME * 0.95) {
                    backtrace_timeout = true;
                    return;
                }

                if (result[i]) {
                    count_one += 1;
                }
            }
            if (!backtrace_timeout && count_one == 0) {
                update_safe(r, c);
            } else if (!backtrace_timeout && count_one == results.size()) {
                update_bomb(r, c);
            } else {
                probabilities.push_back(make_pair((double)count_one/(double)results.size(), p));
            }
        }
    }

    vector<pair<double, Pos>> global_search() {
        if (backtrace_timeout) {
            return vector<pair<double, Pos>>();
        }

        const auto constraints = extract_constraints();
        const auto constraint_groups = split_constraints(constraints);
        vector<pair<double, Pos>> probabilities;
        for (const auto& group : constraint_groups) {
            backtrace(group, probabilities);
            if (backtrace_timeout) {
                dbg(constraint_groups.size());
                break;
            }
        }
        dbg(next_commands.size());
        sort(probabilities.begin(), probabilities.end());
        return probabilities;
    }

    void insert_all_opens() {
        REP(r, params.N) REP(c, params.N) {
            if (grid[r][c].is_hidden() && judge_grid[r][c].is_hidden()) {
                next_commands.push_back(Command::open(r, c));
            }
        }
    }

    void print_stats(string reason) {
        cerr << "reason: " << reason << endl;
        cerr << "runtime: " << get_runtime(true) << endl;
        cerr << "uncover_ratio: " << calc_uncover_ratio() << endl;
        cerr << "mine_hit: " << score_mine_hit << endl;
        stats.print();
        cerr.flush();
    }

    Command decide_next_command() {
        assert(!is_game_end());

        vector<pair<double, Pos>> probabilities;
        if (next_commands.empty()) {
            probabilities = global_search();
        }

        if (grid_bomb_count == params.M && grid_unknown_count > 0 && next_commands.empty()) {
            insert_all_opens();
        }

        if (!next_commands.empty()) {
            Command next = next_commands.back();
            next_commands.pop_back();

            if (next.is_open() && is_last_move()) {
                print_stats("win");
            }
            return next;
        }

        if (backtrace_timeout) {
            print_stats("timeout");
            return Command::stop();
        }

        const int rest = params.M - grid_bomb_count;
        const int all = grid_unknown_count;
        const double default_prob = (double)rest/(double)all;

        function<Command()> factory = [this]() { return this->choose_random_unknown(); };
        double prob = default_prob;
        bool random_guess = true;

        if (!probabilities.empty()) {
            const auto best = probabilities[0];
            dbg(best);
            if (best.first < default_prob) {
                factory = [best]() { return Command::open(best.second.first, best.second.second); };
                prob = best.first;
                random_guess = false;
            }
        }

        const double ratio = calc_uncover_ratio();
        const bool do_invest = score() <= 1.0 / (1.0 + score_mine_hit + ESTIMATED_ADD_BOMB * prob);
        if (do_invest) {

            stats.guess_count += 1;
            if (random_guess) {
                stats.random_guess_count += 1;
            }

            if (is_last_move()) {
                print_stats("maybe_win");
            }

            return factory();
        }

        dbg(score());
        dbg(ratio);
        dbg(score_mine_hit);

        print_stats("stop");
        return Command::stop();
    }

    void fire_update(int row, int col) {
        assert(grid[row][col].is_value());

        const int value = grid[row][col].get_value();
        const int safe_value = neighbors[row][col].size() - value;
        assert(count_safe_neighbors[row][col] <= safe_value);
        assert(count_bomb_neighbors[row][col] <= value);

        if (unknown_neighbors[row][col].empty()) {
            // do nothing.
        } else if (count_safe_neighbors[row][col] == safe_value) {
            const vector<pair<int, int>> nps = vector<pair<int, int>>(unknown_neighbors[row][col].begin(), unknown_neighbors[row][col].end());
            for (auto np : nps) {
                int nr, nc;
                tie(nr, nc) = np;
                if (grid[nr][nc].is_hidden()) {
                    update_bomb(nr, nc);
                }
            }
            positions_with_unknown_values.erase(make_pair(row, col));
        } else if (count_bomb_neighbors[row][col] == value) {
            const vector<pair<int, int>> nps = vector<pair<int, int>>(unknown_neighbors[row][col].begin(), unknown_neighbors[row][col].end());
            for (auto np : nps) {
                int nr, nc;
                tie(nr, nc) = np;
                if (grid[nr][nc].is_hidden()) {
                    update_safe(nr, nc);
                }
            }
            positions_with_unknown_values.erase(make_pair(row, col));
        }
    }

    void update_safe(int row, int col) {
        assert(grid[row][col].is_hidden());
        grid[row][col].set_safe();
        grid_unknown_count -= 1;

        if (judge_grid[row][col].is_hidden()) {
            next_commands.push_back(Command::open(row, col));
        }

        for (auto np : neighbors[row][col]) {
            int nr, nc;
            tie(nr, nc) = np;
            count_safe_neighbors[nr][nc] += 1;
            assert(unknown_neighbors[nr][nc].erase(make_pair(row, col)) == 1);

            if (grid[nr][nc].is_value()) {
                fire_update(nr, nc);
            }
        }
    }

    void update_bomb(int row, int col) {
        assert(grid[row][col].is_hidden());
        grid[row][col].set_bomb();
        grid_unknown_count -= 1;
        grid_bomb_count += 1;

        if (judge_grid[row][col].is_hidden()) {
            next_commands.push_back(Command::flag(row, col));
        }

        for (auto np : neighbors[row][col]) {
            int nr, nc;
            tie(nr, nc) = np;
            count_bomb_neighbors[nr][nc] += 1;
            assert(unknown_neighbors[nr][nc].erase(make_pair(row, col)) == 1);

            if (grid[nr][nc].is_value()) {
                fire_update(nr, nc);
            }
        }
    }

    void update_value(int row, int col, int value) {
        if (grid[row][col].is_hidden()) {
            update_safe(row, col);
        }

        assert(grid[row][col].is_safe());
        grid[row][col].set_value(value);
        if (!unknown_neighbors[row][col].empty()) {
            positions_with_unknown_values.insert(make_pair(row, col));
        }
        // no "grid_unknown_count -= 1" because it's updated in update_safe

        fire_update(row, col);
    }

    void judge_update_value(int row, int col, int value, long runtime) {
        dbg(runtime);
        assert(judge_grid[row][col].is_hidden());
        judge_grid[row][col].set_value(value);
        this->score_uncover += 1;
        set_runtime(runtime);

        if (grid[row][col].is_hidden() || grid[row][col].is_safe()) {
            update_value(row, col, value);
        } else {
            assert(grid[row][col].has_value(value));
        }
        //dbg(score());
    }

    void judge_update_bomb(int row, int col, long runtime) {
        dbg(runtime);
        assert(judge_grid[row][col].is_hidden());
        judge_grid[row][col].set_bomb();
        this->score_mine_hit += 1;
        set_runtime(runtime);

        if (grid[row][col].is_hidden()) {
            update_bomb(row, col);
        } else {
            assert(grid[row][col].is_bomb());
        }
        //dbg(score());
    }
};

int main() {
    string feedback;
    int N, M, D, row, col;
    cin >> N >> M >> D >> row >> col;
    getline(cin, feedback); // read endline

    Solver solver(N, M, D, row, col);

    while(!solver.is_game_end()) {
        Command command = solver.decide_next_command();

        cout << command.to_string() << endl;
        cout.flush();

        if (command.is_stop()) {
            break;
        } else if (command.is_flag()) {
            getline(cin, feedback); // empty line expected
            assert(feedback.empty());
        } else {
            getline(cin, feedback);
            assert(!feedback.empty());
            if (feedback.find("BOOM!")!=string::npos) {
                string dummy;
                long runtime;

                stringstream ss(feedback);
                ss >> dummy >> runtime;
                assert(dummy == "BOOM!");
                solver.judge_update_bomb(command.row(), command.col(), runtime);
            } else {
                int value;
                long runtime;

                stringstream ss(feedback);
                ss >> value >> runtime;
                solver.judge_update_value(command.row(), command.col(), value, runtime);
            }
        }
    }

    return 0;
}
