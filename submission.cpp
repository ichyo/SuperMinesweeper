#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <cassert>

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

    bool has_zero() {
        return value == 0;
    }

    bool is_value() const {
        return value >= 0;
    }

    int get_value() const {
        assert(is_value());
        return value;
    }

    void set_value(int value) {
        assert(is_hidden());
        assert(value >= 0);
        this->value = value;
    }

    void set_bomb() {
        assert(is_hidden());
        this->value = BOMB;
    }
};

class Command {
    enum Type {
        Stop,
        Open,
    };

    Type type;
    int row;
    int col;

    Command(Type type) {
        this->type = type;
    }

    Command(Type type, int row, int col) {
        this->type = type;
        this->row = row;
        this->col = col;
    }

    public:
    static Command stop() {
        return Command(Type::Stop);
    }

    static Command open(int r, int c) {
        return Command(Type::Open, r, c);
    }

    bool is_stop() {
        return type == Type::Stop;
    }

    bool is_open() {
        return type == Type::Open;
    }

    string to_string() {
        if (is_stop()) {
            return "STOP";
        }
        if (is_open()) {
            return "G " + ::to_string(row) + " " + ::to_string(col);
        }
        assert(false);
    }
};

class Solver {
    Params params;

    vector<vector<Cell>> grid;
    int lastR;
    int lastC;
    long runtime;

    bool is_safe[MAX_N][MAX_N];

    bool with_in(int r1, int c1, int r2, int c2, int D) {
        return (r1 - r2) * (r1 - r2) + (c1 - c2) * (c1 - c2) <= D;
    }

    Command prepare_open(int r, int c) {
        this->lastR = r;
        this->lastC = c;
        return Command::open(r, c);
    }

    public:
    Solver(int N, int M, int D, int row, int col) {
        params.N = N;
        params.M = M;
        params.D = D;

        grid = vector<vector<Cell>>(N, vector<Cell>(N, Cell()));
        lastR = row;
        lastC = col;

        REP(r, N) {
            REP(c, N) {
                is_safe[r][c] = false;
            }
        }

        judge_update_value(0, 0);
    }

    Command decide_next_command() {
        REP(r, params.N) {
            REP(c, params.N) {
                if (is_safe[r][c] && grid[r][c].is_hidden()) {
                    return prepare_open(r, c);
                }
            }
        }
        return Command::stop();
    }

    void judge_update_value(int value, long runtime) {
        grid[lastR][lastC].set_value(value);
        this->runtime = runtime;

        if (value == 0) {
            REP(r, params.N) {
                REP(c, params.N) {
                    if (with_in(r, c, lastR, lastC, params.D)) {
                        is_safe[r][c] = true;
                    }
                }
            }
        }
    }

    void judge_update_bomb(long runtime) {
        assert(!is_safe[lastR][lastC]);
        grid[lastR][lastC].set_bomb();
        this->runtime = runtime;
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
        }

        getline(cin, feedback);
        if (feedback.find("BOOM!")!=string::npos) {
            string dummy;
            long runtime;

            stringstream ss(feedback);
            ss >> dummy >> runtime;
            assert(dummy == "BOOM!");
            solver.judge_update_bomb(runtime);
        } else {
            int value;
            long runtime;

            stringstream ss(feedback);
            ss >> value >> runtime;
            solver.judge_update_value(value, runtime);
        }
    }

    return 0;
}
