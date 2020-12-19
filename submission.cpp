#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <string>
#include <cassert>
#include <cstring>
#include <tuple>
#include <algorithm>

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

const int BACKTRACE_MAX_N = 30;

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

class Solver {
    Params params;

    int grid_unknown_count = 0;
    int grid_bomb_count = 0;
    int score_uncover = 0;
    int score_mine_hit = 0;

    vector<vector<Cell>> judge_grid;
    vector<vector<Cell>> grid;

    long runtime = 0;
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

    Command choose_minimum_prob() {
        vector<Pos> probs[MAX_N][MAX_N];
        for(auto p : positions_with_unknown_values) {
            int r, c;
            tie(r, c) = p;
            assert(grid[r][c].is_value());
            assert(!unknown_neighbors[r][c].empty());
            const int total = unknown_neighbors[r][c].size();
            const int rest = grid[r][c].get_value() - count_bomb_neighbors[r][c];
            for(auto np : unknown_neighbors[r][c]) {
                int nr, nc;
                tie(nr, nc) = np;
                probs[nr][nc].push_back(make_pair(rest, total));
            }
        }
        int rest = params.M - grid_bomb_count;
        int all = grid_unknown_count;
        double default_prob = (double)rest/(double)all;

        int best_r = -1, best_c = -1;
        double best_prob = default_prob;

        REP(r, params.N) REP(c, params.N) {
            if (probs[r][c].empty()) {
                continue;
            }
            assert(!probs[r][c].empty());
            double sum = 0;
            for(auto p : probs[r][c]) {
                int a, b;
                tie(a, b) = p;
                sum += (double)a/(double)b;
            }
            const double prob = sum / (double)probs[r][c].size();
            if (prob <= best_prob) {
                best_r = r;
                best_c = c;
                best_prob = prob;
            }
        }
        dbg(best_r, best_c, best_prob);
        if (best_r != -1) {
            return Command::open(best_r, best_c);
        }
        return Command::stop();
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
        vector<bool>& is_mine,
        vector<int>& count_zero,
        vector<int>& count_one,
        vector<vector<bool>>& results,
        const Constraints& constraints,
        const vector<vector<int>>& p2c
    ) {
        if (index == is_mine.size()) {
            results.push_back(is_mine);
            return;
        }
        {
            for(int i : p2c[index]) {
                count_zero[i] ++;
            }
            bool violate = false;
            for(int i : p2c[index]) {
                if (count_zero[i] > constraints[i].max_zero()) {
                    violate = true;
                    break;
                }
            }
            is_mine[index] = false;
            if (!violate) {
                dfs(index + 1, is_mine, count_zero, count_one, results, constraints, p2c);
            }
            for(int i : p2c[index]) {
                count_zero[i] --;
            }
        }
        {
            for(int i : p2c[index]) {
                count_one[i] ++;
            }
            bool violate = false;
            for(int i : p2c[index]) {
                if (count_one[i] > constraints[i].max_one()) {
                    violate = true;
                    break;
                }
            }
            is_mine[index] = true;
            if (!violate) {
                dfs(index + 1, is_mine, count_zero, count_one, results, constraints, p2c);
            }
            for(int i : p2c[index]) {
                count_one[i] --;
            }
        }
    }

    void backtrace(const Constraints& constraints) {
        set<Pos> point_set;
        for(const auto& c : constraints) {
            for(const auto& p : c.positions) {
                point_set.insert(p);
            }
        }
        vector<Pos> points(point_set.begin(), point_set.end());
        if (points.size() > BACKTRACE_MAX_N) {
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
        vector<bool> is_mine(points.size(), false);
        vector<int> count_zero(constraints.size(), 0);
        vector<int> count_one(constraints.size(), 0);

        vector<vector<bool>> results;
        dfs(0, is_mine, count_zero, count_one, results, constraints, p2c);

        assert(!results.empty());

        for (int i = 0; i < points.size(); i++) {
            const auto& p = points[i];
            int r, c;
            tie(r, c) = p;

            if (!grid[r][c].is_hidden()) {
                continue;
            }

            bool first = results[0][i];
            bool all_common = true;
            for(const auto& result : results) {
                if (result[i] != first) {
                    all_common = false;
                    break;
                }
            }
            if (all_common) {
                if (first) {
                    update_bomb(r, c);
                } else {
                    update_safe(r, c);
                }
            }
        }
    }

    void global_search() {
        const auto constraints = extract_constraints();
        const auto constraint_groups = split_constraints(constraints);
        for (const auto& group : constraint_groups) {
            backtrace(group);
        }
    }

    void insert_all_opens() {
        REP(r, params.N) REP(c, params.N) {
            if (grid[r][c].is_hidden() && judge_grid[r][c].is_hidden()) {
                next_commands.push_back(Command::open(r, c));
            }
        }
    }

    Command decide_next_command() {
        if (runtime < MAX_RUNTIME * 0.8 && next_commands.empty()) {
            global_search();
        }

        if (grid_bomb_count == params.M && grid_unknown_count > 0 && next_commands.empty()) {
            insert_all_opens();
        }

        if (!next_commands.empty()) {
            Command next = next_commands.back();
            next_commands.pop_back();
            return next;
        }

        if (calc_uncover_ratio() / (1.0 + score_mine_hit) <= 1.0 / (1.0 + score_mine_hit + 1)) {
            Command next = choose_minimum_prob();
            if (!next.is_stop()) {
                return next;
            }
        }

        if (calc_uncover_ratio() / (1.0 + score_mine_hit) <= 1.0 / (1.0 + score_mine_hit + 1)) {
            Command next = choose_random_unknown();
            if (!next.is_stop()) {
                return next;
            }
        }

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
        this->runtime = runtime;
        this->score_uncover += 1;

        if (grid[row][col].is_hidden() || grid[row][col].is_safe()) {
            update_value(row, col, value);
        } else {
            assert(grid[row][col].has_value(value));
        }
    }

    void judge_update_bomb(int row, int col, long runtime) {
        dbg(runtime);
        assert(judge_grid[row][col].is_hidden());
        judge_grid[row][col].set_bomb();
        this->runtime = runtime;
        this->score_mine_hit += 1;

        if (grid[row][col].is_hidden()) {
            update_bomb(row, col);
        } else {
            assert(grid[row][col].is_bomb());
        }
    }
};

int main() {
    string feedback;
    int N, M, D, row, col;
    cin >> N >> M >> D >> row >> col;
    getline(cin, feedback); // read endline

    Solver solver(N, M, D, row, col);

    while(true) {
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
