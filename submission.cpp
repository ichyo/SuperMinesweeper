#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <cassert>
#include <cstring>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <chrono>
#include <bitset>
#include <functional>
#include <stdint.h>

#ifdef LOCAL_TEST
#include "./dbg.hpp"
#define MAX_RUNTIME 7500
#else
#define dbg(...)
#define MAX_RUNTIME 10000
#endif

#define REP(i,n) for(int i = 0; i < (n); ++i)

using namespace std;

/*  Written in 2019 by David Blackman and Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */


/* This is xoshiro256++ 1.0, one of our all-purpose, rock-solid generators.
   It has excellent (sub-ns) speed, a state (256 bits) that is large
   enough for any parallel application, and it passes all tests we are
   aware of.

   For generating just floating-point numbers, xoshiro256+ is even faster.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s. */

namespace rng {

    /*  Written in 2015 by Sebastiano Vigna (vigna@acm.org)

        To the extent possible under law, the author has dedicated all copyright
        and related and neighboring rights to this software to the public domain
        worldwide. This software is distributed without any warranty.

        See <http://creativecommons.org/publicdomain/zero/1.0/>. */


    /* This is a fixed-increment version of Java 8's SplittableRandom generator
       See http://dx.doi.org/10.1145/2714064.2660195 and
       http://docs.oracle.com/javase/8/docs/api/java/util/SplittableRandom.html

       It is a very fast generator passing BigCrush, and it can be useful if
       for some reason you absolutely want 64 bits of state. */

    namespace splitmix64 {
        static uint64_t x; /* The state can be seeded with any value. */

        uint64_t next() {
            uint64_t z = (x += 0x9e3779b97f4a7c15);
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
            z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
            return z ^ (z >> 31);
        }
    }
    static inline uint64_t rotl(const uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }

    static uint64_t s[4];

    void seed(uint64_t seed) {
        splitmix64::x = seed;
        s[0] = splitmix64::next();
        s[1] = splitmix64::next();
        s[2] = splitmix64::next();
        s[3] = splitmix64::next();
    }

    uint64_t next(void) {
        const uint64_t result = rotl(s[0] + s[3], 23) + s[0];

        const uint64_t t = s[1] << 17;

        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];

        s[2] ^= t;

        s[3] = rotl(s[3], 45);

        return result;
    }

    /* This is the jump function for the generator. It is equivalent
       to 2^128 calls to next(); it can be used to generate 2^128
       non-overlapping subsequences for parallel computations. */

    void jump(void) {
        static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };

        uint64_t s0 = 0;
        uint64_t s1 = 0;
        uint64_t s2 = 0;
        uint64_t s3 = 0;
        for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
            for(int b = 0; b < 64; b++) {
                if (JUMP[i] & UINT64_C(1) << b) {
                    s0 ^= s[0];
                    s1 ^= s[1];
                    s2 ^= s[2];
                    s3 ^= s[3];
                }
                next();
            }

        s[0] = s0;
        s[1] = s1;
        s[2] = s2;
        s[3] = s3;
    }



    /* This is the long-jump function for the generator. It is equivalent to
       2^192 calls to next(); it can be used to generate 2^64 starting points,
       from each of which jump() will generate 2^64 non-overlapping
       subsequences for parallel distributed computations. */

    void long_jump(void) {
        static const uint64_t LONG_JUMP[] = { 0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe635 };

        uint64_t s0 = 0;
        uint64_t s1 = 0;
        uint64_t s2 = 0;
        uint64_t s3 = 0;
        for(int i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
            for(int b = 0; b < 64; b++) {
                if (LONG_JUMP[i] & UINT64_C(1) << b) {
                    s0 ^= s[0];
                    s1 ^= s[1];
                    s2 ^= s[2];
                    s3 ^= s[3];
                }
                next();
            }
        s[0] = s0;
        s[1] = s1;
        s[2] = s2;
        s[3] = s3;
    }
}

const int MAX_N = 50;

const int CANCEL_LIMIT = MAX_RUNTIME / 14;

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
    double max_score = 0.0;

    void print() {
        cerr << "guess_count: " << guess_count << endl;
        cerr << "random_guess_count: " << random_guess_count << endl;
        cerr << "max_score: " << max_score << endl;
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
            int r = rng::next() % params.N;
            int c = rng::next() % params.N;
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

    vector<pair<int, Pos>> find_new_neighbors() {
        set<Pos> all_unknown_neighbors;
        for(auto p : positions_with_unknown_values) {
            int r, c;
            tie(r, c) = p;
            assert(grid[r][c].is_value() && !unknown_neighbors[r][c].empty());
            for(auto np : unknown_neighbors[r][c]) {
                all_unknown_neighbors.insert(np);
            }
        }
        map<Pos, int> counter;
        for(auto p : all_unknown_neighbors) {
            int r, c;
            tie(r, c) = p;
            for(auto np : unknown_neighbors[r][c]) {
                if (!all_unknown_neighbors.count(np)) {
                    counter[np] += 1;
                }
            }
        }
        vector<pair<int, Pos>> result;
        for(auto p : counter) {
            result.push_back(make_pair(p.second, p.first));
        }
        sort(result.begin(), result.end());
        reverse(result.begin(), result.end());
        return result;
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

    vector<char> _dfs_is_mine;
    vector<int> _dfs_count_zero;
    vector<int> _dfs_count_one;
    vector<int> _dfs_result_ones;
    int _dfs_result_total;
    Constraints _dfs_constraints;
    vector<vector<int>> _dfs_p2c;
    int _dfs_first_one_index;
    bool _dfs_cancel;

    void dfs(int index){
        if (index == _dfs_p2c.size()) {
            _dfs_result_total += 1;
            for (int i = 0; i < _dfs_result_ones.size(); i++) {
                if (_dfs_is_mine[i]) {
                    _dfs_result_ones[i] += 1;
                }
            }
            return;
        }

        if (backtrace_timeout) {
            return;
        }

        if (_dfs_cancel) {
            return;
        }

        if (get_runtime() > MAX_RUNTIME * 0.95) {
            dbg(get_runtime());
            dbg(_dfs_p2c.size());
            dbg(_dfs_constraints.size());
            dbg(_dfs_result_total);
            backtrace_timeout = true;
            return;
        }
        if (params.D >= 8 && get_runtime() - runtime > CANCEL_LIMIT) {
            _dfs_cancel = true;
            return;
        }

        {
            int violate = _dfs_p2c[index].size();
            for (int i = 0; i < _dfs_p2c[index].size(); i++) {
                const int x = _dfs_p2c[index][i];
                if (++_dfs_count_zero[x] > _dfs_constraints[x].max_zero()) {
                    violate = i;
                    break;
                }
            }
            if (violate == _dfs_p2c[index].size()) {
                _dfs_is_mine[index] = false;
                dfs(index + 1);
                violate--;
            }
            for (int i = 0; i <= violate; i++) {
                const int x = _dfs_p2c[index][i];
                --_dfs_count_zero[x];
            }
        }
        {
            if (index < _dfs_first_one_index) {
                _dfs_first_one_index = index;
            }

            int violate = _dfs_p2c[index].size();
            for (int i = 0; i < _dfs_p2c[index].size(); i++) {
                const int x = _dfs_p2c[index][i];
                if (++_dfs_count_one[x] > _dfs_constraints[x].max_one()) {
                    violate = i;
                    break;
                }
            }
            if (violate == _dfs_p2c[index].size()) {
                _dfs_is_mine[index] = true;
                dfs(index + 1);
                violate--;
            }
            for (int i = 0; i <= violate; i++) {
                const int x = _dfs_p2c[index][i];
                --_dfs_count_one[x];
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

        map<Pos, int> occurance;
        for(const auto& c : constraints) {
            for(const auto& p : c.positions) {
                occurance[p]++;
            }
        }

        sort(points.begin(), points.end(), [&](const Pos& a, const Pos& b){
            return occurance[a] > occurance[b];
        });
        assert(occurance[points[0]] >= occurance[points[points.size() - 1]]);

        _dfs_p2c = vector<vector<int>>(points.size(), vector<int>());
        for(int i = 0; i < constraints.size(); i++) {
            for (int j = 0; j < points.size(); j++) {
                if (constraints[i].positions.count(points[j])) {
                    _dfs_p2c[j].push_back(i);
                }
            }
        }

        this->_dfs_constraints = constraints;
        _dfs_is_mine = vector<char>(points.size());
        _dfs_count_zero = vector<int>(constraints.size(), 0);
        _dfs_count_one = vector<int>(constraints.size(), 0);
        _dfs_result_ones = vector<int>(points.size());
        _dfs_result_total = 0;
        _dfs_first_one_index = points.size();
        _dfs_cancel = false;
        dfs(0);

        if (backtrace_timeout) {
            return;
        }

        if (_dfs_cancel) {
            probabilities.push_back(make_pair(0.2, points[0])); // make to pick up first
            return;
        }

        assert(_dfs_result_total > 0);

        for (int i = 0; i < points.size(); i++) {
            const auto& p = points[i];
            int r, c;
            tie(r, c) = p;

            if (!grid[r][c].is_hidden()) {
                continue;
            }

            const int count_one = _dfs_result_ones[i];
            const int count_all = _dfs_result_total;
            if (!backtrace_timeout && count_one == 0) {
                update_safe(r, c);
            } else if (!backtrace_timeout && count_one == count_all) {
                update_bomb(r, c);
            } else {
                probabilities.push_back(make_pair((double)count_one/(double)count_all, p));
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
            dbg(_dfs_result_total);
            if (backtrace_timeout) {
                dbg(constraint_groups.size());
                break;
            }
        }
        dbg(next_commands.size());
        //dbg(probabilities.size());
        dbg(get_runtime(true) - runtime);
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

    double find_maximum_score(
        const int all,
        const int rest,
        const int M, // mine_hit
        const double g
    ){
        const int cur = all - rest;
        const auto f = [&](double x){
            double s1 = (cur + x) / all;
            double s2 = 1.0 / (1.0 + M + g * x);
            return s1 * s2;
        };
        double low = 0;
        double high = rest;
        REP(_, 100) {
            const double c1 = (low * 2 + high) / 3;
            const double c2 = (low + high * 2) / 3;
            if (f(c1) < f(c2)) {
                low = c1;
            } else {
                high = c2;
            }
        }
        return f(low);
    }

    double estimate_best_score(double prob) {
        map<int, tuple<double, double, double>> param_map;
        param_map[1] = make_tuple(0.24931335, 17.86378303, -2.76674817);
        //param_map[2] = make_tuple(4.52466102e-03, 2.69382644e+01, -1.37446064e-01);

        map<int, double> k_map;
        k_map[1] = 0.5;
        //k_map[2] = 0.5;

        if (param_map.count(params.D)) {
            const int all = params.N * params.N - params.M;
            const int cur = score_uncover;
            const int rest = all - cur;

            const double r = (double)params.M/(double)params.N/(double)params.N; // ratio (0.1~0.3)

            const double a = get<0>(param_map[params.D]);
            const double b = get<1>(param_map[params.D]);
            const double c = get<2>(param_map[params.D]);

            const double g = (a * exp(b * r) + c) / 1000.0;

            const double k = k_map[params.D];

            const double maximum_score_good = 1.0 / (1.0 + score_mine_hit + g * rest * k);
            const double maximum_score_bad = 1.0 / (1.0 + score_mine_hit + 1 + g * rest * k);

            //const double maximum_score_good = find_maximum_score(all, rest, score_mine_hit, g * k);
            //const double maximum_score_bad = find_maximum_score(all, rest, score_mine_hit + 1, g * k);

            const double result = prob * maximum_score_bad + (1 - prob) * maximum_score_good;

            dbg(score());
            dbg(maximum_score_good);
            dbg(maximum_score_bad);
            dbg(result);

            return result;
        }

        return 1.0 / (1.0 + score_mine_hit + ESTIMATED_ADD_BOMB * prob);
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
                this->score_uncover += 1; // to print correct stats
                stats.max_score = max(score(), stats.max_score);
                print_stats("win");
                this->score_uncover -= 1;
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

        function<Command()> factory = [this]() {
            return this->choose_random_unknown();
        };
        double prob = default_prob;
        bool random_guess = true;

        if (!probabilities.empty()) {
            const auto best = probabilities[0];
            dbg(best);
            dbg(default_prob);
            if (best.first <= default_prob) {
                factory = [best]() { return Command::open(best.second.first, best.second.second); };
                prob = best.first;
                random_guess = false;
            }
            dbg(random_guess);
        }

        const double ratio = calc_uncover_ratio();
        const bool do_invest = score() <= estimate_best_score(prob);
        if (do_invest) {

            stats.guess_count += 1;
            if (random_guess) {
                stats.random_guess_count += 1;
            }

            if (is_last_move()) {
                this->score_uncover += 1; // to print correct stats
                stats.max_score = max(score(), stats.max_score);
                print_stats("maybe_win");
                this->score_uncover -= 1;
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
        //dbg(runtime);
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

        stats.max_score = max(score(), stats.max_score);
    }

    void judge_update_bomb(int row, int col, long runtime) {
        cerr << "BOMB!" << endl;
        //dbg(runtime);
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

//#include <sys/random.h>
int main() {
    //uint64_t seed;
    //getrandom(&seed, 8, 0);

    rng::seed(2020);

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
