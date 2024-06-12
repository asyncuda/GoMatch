#include <iostream>
#include <cmath>
#include <algorithm>
#include <list>
#include <utility>
#include <vector>
#include <climits>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <cmath>

struct member { int id; std::string name; int rank; std::vector<int> score; };

using mpair = std::pair<member, member>;

using group = std::vector<mpair>;

std::vector<int> exclusion_table;

std::vector<std::pair<int,int>> matched_table;

const int w1 = 1;

const int w2 = 1000;

const int w3 = 100;

const int w4 = 200;

int N = 0;

int diff_rank(const mpair &mp) { return std::abs(mp.first.rank - mp.second.rank); }

int played(const mpair &mp) { return matched_table[mp.first.id * N + mp.second.id].first + matched_table[mp.first.id * N + mp.second.id].second; }

double win_rate(int id) {
    int win = 0, lose = 0;
    for(auto itr = matched_table.begin() + id * N; itr != std::next(matched_table.begin() + id * N, N); ++itr) {
        win += itr->first;
        lose += itr->second;
    }
    if(win + lose == 0) {
        return 0.;
    }

    return (double)win/(win+lose);
}

int win_count(int id) 
{
    int sum = 0;
    for(int i = 0; i < N; i++) { sum += matched_table[id * N + i].first; }
    return sum;
}

int sos(int id)
{
    int sum = 0;
    for(int i = 0 ; i < N; i++) {
        if(matched_table[id * N + i].first + matched_table[id * N + i].second > 0) {
            sum += win_count(i);
        }
    }
    return sum;
}

int sosos(int id)
{
    int sum = 0;
    for(int i = 0 ; i < N; i++) {
        if(matched_table[id * N + i].first + matched_table[id * N + i].second > 0) {
            sum += sos(i);
        }
    }
    return sum;
}

double diff_score(const mpair &mp) { return std::abs(win_rate(mp.first.id) - win_rate(mp.second.id)); }

int exclusion(const mpair &mp) { return exclusion_table[mp.first.id * N + mp.second.id]; }

//  ∑ (w1*棋力差)^2 + (w2*対局済)^2 + (w3*勝率の差)^2 + (w4*避けたい組合せ)^2)
int hamiltonian(group &match)
{
    return std::accumulate(match.begin(), match.end(), 0, [](int acc, mpair mp) {
            return acc + std::pow(w1 * diff_rank(mp), 2.) + std::pow(w2 * played(mp), 2.) + std::pow(w3 * diff_score(mp), 2.) + std::pow(w4 * exclusion(mp), 2.); }
        );
}

group dfs(std::list<member> rest, group candidate);

void get_members(std::list<member> *members, std::string filename);

void import_exclusion(std::string filename);

void import_table(std::string filename);

void export_table(std::string filename, group &match);

int main(int argc, char *argv[])
{
    if(argc != 2) {
        std::cerr << "useage: ./a.out block-name" << std::endl;
        std::exit(-1);
    }

    std::string blockname(argv[1]);

    std::list<member> members;

    get_members(&members, blockname);

    matched_table.resize(N*N);
    exclusion_table.resize(N*N);
    for(int i = 0; i < N*N; i++) {
        matched_table[i] = std::make_pair(0, 0);
        exclusion_table[i] = 0;
    }

    import_table(blockname);
    import_exclusion(blockname);

    std::cout << "id name rank win_rate S SOS" << std::endl;
    for(auto &x : members) {
        x.score.push_back(win_count(x.id));
        x.score.push_back(sos(x.id));
        x.score.push_back(sosos(x.id));
    }

    members.sort([](member &lhs, member &rhs) {
            for(int i = 0; i < lhs.score.size(); i++) {
                if(lhs.score[i] != rhs.score[i]) {
                    return lhs.score[i] > rhs.score[i];
                }
            }
            return lhs.score[0] > rhs.score[0];
            });

    for(auto &x : members) {
        std::cout << x.id << "\t" << x.name << "\t" << x.rank << "\t" << win_rate(x.id) << "\t" << x.score[0] << "\t" << x.score[1] << "\t" << x.score[2] << std::endl;
    }

    auto print_match = [&](const group gr) {
        for(auto g : gr) { std::cout << g.first.name << "(" << g.first.id << ")vs" << g.second.name << "(" << g.second.id << ")" << std::endl; 
            if(played(g)) {
                std::cout << "already!" << played(g) << " ";

            }
        } 
        std::cout << std::endl;
    };

    group best;

    std::cout << std::setprecision(2);

    best = dfs(members, group());
    print_match(best);

    export_table(blockname, best);

    return 0;
}

group dfs(std::list<member> rest, group candidate) 
{
    if(rest.empty()) {
        return candidate;
    }

    member f = rest.front(); rest.pop_front();
    std::pair<int, group> best; best.first = INT_MAX;
    for(auto itr = rest.begin(); itr != rest.end(); ++itr) {
        member s = *itr;

        auto pos = rest.erase(itr);
        candidate.push_back(std::make_pair(f, s));

        auto tmp = dfs(rest, candidate);
        if(hamiltonian(tmp) < best.first) {
            best = make_pair(hamiltonian(tmp), tmp);
        }

        candidate.pop_back();
        itr = rest.insert(pos, s);
    }

    return best.second;
}

void get_members(std::list<member> *members, std::string filename)
{
    std::string name;
    int rank = 0;

    std::ifstream ifs(filename + ".member");

    ifs >> N;

    for(int i = 0; i < N; i++) {
        ifs >> name >> rank;
        members->push_back({i, name, rank});
    }
}

void import_exclusion(std::string filename)
{
    std::ifstream ifs(filename + ".exclusion");

    int len; ifs >> len;

    int x, y;
    for(int i = 0; i < len; i++) {
        ifs >> x >> y;
        exclusion_table[x*N+y] = 1;
        exclusion_table[y*N+x] = 1;
    }
}

void import_table(std::string filename)
{
    std::ifstream ifs(filename + ".table");

    int len; ifs >> len;

    int myid, opid, win, lose;

    for(int i = 0; i < len; i++) {
        ifs >> myid >> opid >> win >> lose;
        matched_table[myid * N + opid].first = win;
        matched_table[myid * N + opid].second = lose;
        matched_table[opid * N + myid].first = lose;
        matched_table[opid * N + myid].second = win;
    }
}

void export_table(std::string filename, group &match)
{
    std::vector<std::string> dat;

    for(auto mp : match) {
        auto left = std::min(mp.first.id, mp.second.id);
        auto right= std::max(mp.first.id, mp.second.id);
        dat.push_back(std::to_string(left) + " " + std::to_string(right) + " " + 
                std::to_string(matched_table[left*N+right].first) + " " + std::to_string(matched_table[left*N+right].second) + " ?");
    }

    auto check = [&](int i, int j) {
        for(auto mp : match) {
            if(i == std::min(mp.first.id, mp.second.id) and j == std::max(mp.first.id, mp.second.id)) {
                return true;
            }
        }
        return false;
    };

    for(int i = 0; i < N; i++) {
        for(int j = i + 1; j < N; j++) {
            if(matched_table[i*N+j].first == 0 and matched_table[i*N+j].second == 0) {
                continue;
            }
            if(check(i, j)) {
                continue;
            }
            dat.push_back(std::to_string(i) + " " + std::to_string(j) + " " + 
                    std::to_string(matched_table[i*N+j].first) + " " + std::to_string(matched_table[i*N+j].second));
        }
    }

    std::ofstream ofs(filename + ".table");
    ofs << dat.size() << std::endl;
    for(int i = 0; i < dat.size(); i++) {
        ofs << dat[i] << std::endl;
    }
}
