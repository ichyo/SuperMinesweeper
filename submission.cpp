#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <string>
#include <cassert>
#include <cstring>
#include <tuple>

#define REP(i,n) for(int i = 0; i < (n); ++i)

using namespace std;

const int MAX_N = 50;

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

    vector<vector<Cell>> judge_grid;
    vector<vector<Cell>> grid;

    long runtime;
    vector<Command> next_commands;

    int count_safe_neighbors[MAX_N][MAX_N];
    int count_bomb_neighbors[MAX_N][MAX_N];
    set<pair<int, int>> unknown_neighbors[MAX_N][MAX_N];
    vector<pair<int, int>> neighbors[MAX_N][MAX_N];

    bool with_in(int r1, int c1, int r2, int c2, int D) {
        return (r1 - r2) * (r1 - r2) + (c1 - c2) * (c1 - c2) <= D;
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

        judge_update_value(row, col, 0, 0);
    }

    Command decide_next_command() {
        if (!next_commands.empty()) {
            Command next = next_commands.back();
            next_commands.pop_back();
            return next;
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
        } else if (count_bomb_neighbors[row][col] == value) {
            const vector<pair<int, int>> nps = vector<pair<int, int>>(unknown_neighbors[row][col].begin(), unknown_neighbors[row][col].end());
            for (auto np : nps) {
                int nr, nc;
                tie(nr, nc) = np;
                if (grid[nr][nc].is_hidden()) {
                    update_safe(nr, nc);
                }
            }
        }
    }

    void update_safe(int row, int col) {
        assert(grid[row][col].is_hidden());
        grid[row][col].set_safe();

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

        fire_update(row, col);
    }

    void judge_update_value(int row, int col, int value, long runtime) {
        assert(judge_grid[row][col].is_hidden());
        judge_grid[row][col].set_value(value);
        this->runtime = runtime;

        if (grid[row][col].is_hidden() || grid[row][col].is_safe()) {
            update_value(row, col, value);
        } else {
            assert(grid[row][col].has_value(value));
        }
    }

    void judge_update_bomb(int row, int col, long runtime) {
        assert(judge_grid[row][col].is_hidden());
        judge_grid[row][col].set_bomb();
        this->runtime = runtime;

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
        } else {
            getline(cin, feedback);
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
