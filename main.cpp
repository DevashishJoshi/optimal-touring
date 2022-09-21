#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <string>
#include <sstream>
#include <limits.h>
#include <time.h>
#include <math.h>

using namespace std;

const int lower_cluster_bound = 3, upper_cluster_bound = 30, repeat = 1, mid_hour_multiplier = 5, time_limit = 60, epsilon = 0.1; 
const bool debug = false;

class City{
    public:
        int street, avenue, mid_hour = -1;
        int id = -1, desired_time, begin_hour = -1, end_hour = -1, value, cluster_number;

        City() { }
        City(int s, int a, int st, int d, int v) : id(s), avenue(a), street(st), desired_time(d), value(v) {}
};

unordered_map<int, unordered_map<int, City>> days_to_city;
unordered_map<int, City> cities_hash_map;
unordered_map<int, vector<pair<int, int>>> cities_day_cluster;
unordered_map<int, unordered_map<int, pair<pair<int, int>, int>>> cluster_centers;
unordered_map<int, int> cluster_curr_time;
unordered_map<int, unordered_map<int, unordered_set<int>>> cluster_unvisited_cities;
unordered_map<int, unordered_map<int, unordered_map<int, City>>> days_to_cluster;
unordered_set<int> visited;

double get_dist(const pair<pair<int, int>, int> &p1, pair<pair<int, int>, int> &p2){
    return abs(p1.first.first - p2.first.first) + abs(p1.first.second - p2.first.second) + abs(p1.second - p2.second);
}

double get_dist(const pair<pair<int, int>, int> &p1, City &city){
    return abs(p1.first.first - city.avenue) + abs(p1.first.second - city.street) + abs(p1.second - city.mid_hour);
}

double get_dist(City &city1, City &city2){
    return abs(city1.avenue - city2.avenue) + abs(city1.street - city2.street) + abs(city1.mid_hour - city2.mid_hour);
}

pair<pair<int, int>, int> get_cluster_center(unordered_map<int, City> &cities){
    double c_x = 0, c_y = 0, c_z = 0;
    int n = 0;

    for(auto &it: cities){
        auto &city = it.second;
        c_x += city.avenue; c_y += city.street; c_z += city.mid_hour;
        n++;
    }

    c_x /= n; c_y /= n; c_z /= n;

    return make_pair(make_pair(c_x, c_y), c_z);
}

bool is_eligible(City &city, int day, int cluster_number){
    bool is_eligible_ = (max(city.begin_hour, cluster_curr_time[day]) + city.desired_time) <= city.end_hour;

    is_eligible_ *= cluster_unvisited_cities[day][cluster_number].find(city.id) != cluster_unvisited_cities[day][cluster_number].end();

    is_eligible_ *= visited.find(city.id) == visited.end();

    return is_eligible_;
}

double calculate_goodness(City &city, int day, int cluster_number){
    return city.value * city.value * 1.0 / (max(city.begin_hour, cluster_curr_time[day]) + city.desired_time);
}

int get_next_city(int day, int cluster_number, int prev_city = -1){
    vector<int> remove_me;

    int best_city = -1, highest_goodness = 0;
    for(int id: cluster_unvisited_cities[day][cluster_number]){
        auto &city_ = days_to_city[day][id];

        if(is_eligible(city_, day, cluster_number)){
            auto goodness = calculate_goodness(city_, day, cluster_number);

            if(prev_city != -1){
                auto &city_prev = days_to_city[day][prev_city];
                goodness *= sqrt(get_dist(city_, city_prev));
            }

            if(goodness > highest_goodness){
                highest_goodness = goodness;
                best_city = id;
            }
        }
        else{
            remove_me.push_back(id);
        }
    }

    for(auto &r: remove_me){
        cluster_unvisited_cities[day][cluster_number].erase(r);
        days_to_city[day].erase(r);
        days_to_cluster[day][cluster_number].erase(r);
    }

    return best_city;
}

int select_cluster(unordered_map<int, unordered_map<int, City>> &clusters, int day, int previous_cluster = -1){
    int best_cluster = -1, best_cluster_goodness = 0;
    pair<pair<int, int>, int> center2;

    if(previous_cluster != -1){
        center2 = cluster_centers[day][previous_cluster];
    }

    for(auto &it: clusters){
        int cluster_number = it.first;
        auto &cluster = it.second;
        
        vector<double> goodness;
        // TODO: Fix this!!
        //vector<int> remove_me;

        for(auto &it2: cluster){
            auto &city = it2.second;
            if(is_eligible(city, day, cluster_number)){
                goodness.push_back(calculate_goodness(city, day, cluster_number));
            }
            else{
                //remove_me.push_back(city.id);
            }
        }

        /*for(auto &r: remove_me){
            cluster_unvisited_cities[day][cluster_number].erase(r);
            days_to_city[day].erase(r);
            days_to_cluster[day][cluster_number].erase(r);
        }*/

        if(goodness.size() == 0){
            // This cluster doesn't have any eligible cities, move on now
            continue;
        }

        auto avg_goodness = accumulate(goodness.begin(), goodness.end(), 0.0) / goodness.size();
        auto center1 = cluster_centers[day][cluster_number];
        int count_of_taken = 0;

        for(auto &it3: days_to_cluster[day][cluster_number]){
            auto &city = it3.second;
            if(cities_day_cluster.find(city.id) == cities_day_cluster.end()){
                count_of_taken++;
            }
        }

        if(previous_cluster != -1){
            avg_goodness *= sqrt(get_dist(center1, center2) / (epsilon + count_of_taken));
        }

        if(best_cluster_goodness < avg_goodness){
            best_cluster_goodness = avg_goodness;
            best_cluster = cluster_number;
        }
    }

    return best_cluster;
}

int get_nearest(unordered_map<int, pair<pair<int, int>, int>> &centers, City &city){
    // TODO: Improve this!!!!!
    if(centers.empty()) return -1;

    double nearest_dist = (double)INT_MAX;
    int nearest = 0;
    for(auto &it: centers){
        double dist = get_dist(it.second, city);
        if(nearest_dist < dist){
            nearest_dist = dist;
            nearest = it.first;
        }
    }

    return nearest;
}

double get_silhouette_score(int day, int n_clusters){
    vector<pair<double, double>> a_and_b_s;
    // a is the mean intra-cluster distance, b is the mean inter-cluster distance

    double a, a_n, avg = 0;
    auto &cities = days_to_city[day];
    for(auto &it : cities){
        auto &city = it.second;
        a = 0; a_n = -1;
        unordered_map<int, pair<double, int>> b;
        for(auto &it2: cities){
            auto &other = it2.second;
            if(other.cluster_number == city.cluster_number){
                a += get_dist(days_to_city[day][other.id], days_to_city[day][city.id]);
                a_n++;
            }
            else{
                if(b.find(city.cluster_number) == b.end()){
                    b[city.cluster_number].first = 0;
                    b[city.cluster_number].second = 0;
                }

                b[city.cluster_number].first += get_dist(days_to_city[day][other.id], days_to_city[day][city.id]);
                b[city.cluster_number].second++;
            }
        }

        a /= a_n;
        double b_final = (double)INT_MAX;
        for(auto &it: b){
            b_final = min(b_final, it.second.first / it.second.second);
        }

        a_and_b_s.push_back(make_pair(a, b_final));
    }
    
    for(auto &it: a_and_b_s){
        a = it.first; double b_final = it.second;
        avg += (b_final - a) / max(a, b_final);
    }

    return avg /= a_and_b_s.size();
}

const unordered_map<int, unordered_map<int, City>> make_clusters(unordered_map<int, City> &cities, int day, int n_clusters){
    unordered_map<int, unordered_map<int, City>> clusters;

    // Assign each city a cluster number from 1 to n_clusters
    for (auto &it: cities){
        auto &city = it.second;
        city.cluster_number = rand() % n_clusters + 1;
        clusters[city.cluster_number][city.id] = city;
    }

    unordered_map<int, pair<pair<int, int>, int>> centers;

    int n;
    double c_x, c_y, c_z;
    bool converged = false;
    while(!converged){
        converged = true;
        for(auto &it: clusters){
            c_x = 0; c_y = 0; c_z = 0;
            n = 0;
            for(auto &it2: it.second){
                auto &city = it2.second;

                c_x += city.avenue; c_y += city.street; c_z += city.mid_hour;
                n++;
            }

            c_x /= n; c_y /= n; c_z /= n;
            centers[it.first] = make_pair(make_pair(c_x, c_y), c_z);
        }

        for(auto &cluster: clusters){
            for(auto &it3: cluster.second){
                auto &city = it3.second;
                int new_cluster = get_nearest(centers, city);

                if(new_cluster != city.cluster_number){
                    converged = false;
                }

                city.cluster_number = new_cluster;
            }   
        }
    }

    return clusters;
}

const unordered_map<int, unordered_map<int, City>> make_clusters(unordered_map<int, City> &cities, int day){
    int n_clusters = lower_cluster_bound, best = n_clusters;
    double max_silhouette_score = -1;
    unordered_map<int, unordered_map<int, City>> return_me;
    int repeat_count = 0;
    
    while(n_clusters < min((int)cities.size(), upper_cluster_bound)){
        auto &clusters = make_clusters(cities, day, n_clusters);
        double silhouette_score = get_silhouette_score(day, n_clusters);

        if(max_silhouette_score < silhouette_score){
            max_silhouette_score = silhouette_score;
            best = n_clusters;
            return_me = clusters;
        }

        repeat_count = (repeat_count + 1) % repeat;
        if(repeat_count == 0){
            n_clusters++;
        }
    }

    return return_me;
}

int main(){
    auto start_time = time(NULL);

    string line;
    ifstream infile("/mnt/c/Users/Dev/HPS/optimal-touring/input.txt");

    getline(infile, line);
    while (getline(infile, line)){
        istringstream iss(line);
        int s, a, st, d, v;
        if (!(iss >> s >> a >> st >> d >> v)) { break; }
        cities_hash_map[s] = City(s, a, st, d, v);
    }

    getline(infile, line);
    while (getline(infile, line)){
        istringstream iss(line);
        int s, d, b, e;
        if (!(iss >> s >> d >> b >> e)) { break; }

        City c = City(cities_hash_map[s]);

        c.begin_hour = b * 60;
        c.end_hour = e * 60;
        c.mid_hour = (b + e) * mid_hour_multiplier / 2;

        days_to_city[d][s] = c;
    }

    auto n_days = days_to_city.size();

    vector<vector<int>> final_ans(n_days + 1);

    if(debug)
        cout << "Done reading inputs" << endl;

    // Time to run k means for each day
    for (auto &it: days_to_city){
        int day = it.first;
        auto &cities = it.second;
        days_to_cluster[it.first] = make_clusters(cities, day);
    }

    if(debug)
        cout << "Done creating clusters" << endl;

    // Clear cluster numbers
    for(auto &it: days_to_city){
        for(auto &it2: it.second){
            it2.second.cluster_number = -1;
        }
    }

    if(debug)
        cout << "Done clearing cluster numbers" << endl;

    // Find final cluster centers
    for(auto &it: days_to_cluster){
        int day = it.first;
        for(auto &it2: it.second){
            int cluster_number = it2.first;
            auto &clusters = it2.second;
            cluster_centers[day][cluster_number] = get_cluster_center(clusters);
            cluster_curr_time[day] = 0;
        }
    }

    if(debug)
        cout << "Found cluster centers" << endl;

    // Fill the hash set of unvisited cities
    for(auto &it: days_to_cluster){
        int day = it.first;
        for(auto &it2: it.second){
            int cluster_number = it2.first;
            auto &clusters = it2.second;
            for(auto &it: clusters){
                auto &city = it.second;

                cluster_unvisited_cities[day][cluster_number].insert(city.id);
                cities_day_cluster[city.id].push_back(make_pair(day, cluster_number));
            }
        }
    }

    if(debug)
        cout << "Filled the hash set of unvisited cities" << endl;

    // Now keep selecting a good cluster, and visit cities in that cluster alternating for each day
    unordered_map<int, int> current_cluster;
    unordered_map<int, int> previous_cluster;

    unordered_map<int, unordered_map<int, int>> previous_city;

    unordered_set<int> done;

    if(debug)
        cout << "Will go through the while loop now" << endl;

    while((time(NULL) - start_time) <= time_limit && done.size() < n_days){
        for(auto &it: days_to_cluster){
            int day = it.first;

            auto &clusters = it.second;
            if(done.find(day) != done.end()) continue;

            int cluster_number;

            if(current_cluster.find(day) == current_cluster.end()){
                // There is no current cluster for this day, will select a good cluster now
                if(previous_cluster.find(day) != previous_cluster.end()){
                    cluster_number = select_cluster(clusters, day, previous_cluster[day]);
                }
                else{
                    cluster_number = select_cluster(clusters, day);
                }
                
                if(cluster_number == -1){
                    // No clusters left, so we move on to the next day
                    done.insert(day);
                    continue;
                }
                current_cluster[day] = cluster_number;
            }
            else{
                cluster_number = current_cluster[day];
            }

            // Find a good valid city in the current cluster for this day. If there is none, then skip a turn
            // and clear current_cluster's entry
            // But if there is one, then visit it and move on to the next day
            
            int good_city = -1;
            if(previous_city.find(day) == previous_city.end() || previous_city[day].find(cluster_number) == previous_city[day].end()){
                good_city = get_next_city(day, cluster_number);
            }
            else{
                good_city = get_next_city(day, cluster_number, previous_city[day][cluster_number]);
            }

            if(debug && good_city != -1)
                cout << "good_city: " << good_city << endl;

            previous_city[day][cluster_number] = good_city;

            if(good_city == -1){
                // Time to move on to the next cluster, this one is all done
                previous_cluster[day] = cluster_number;
                current_cluster.erase(day);
                continue;
            }
            else{
                for(auto &day_cluster: cities_day_cluster[good_city]){
                    int day = day_cluster.first;
                    int cluster = day_cluster.second;

                    cluster_unvisited_cities[day][cluster].erase(good_city);
                }
                
                final_ans[day].push_back(good_city);
                visited.insert(good_city);
                cluster_curr_time[day] = max(cities_hash_map[good_city].begin_hour, cluster_curr_time[day]) + cities_hash_map[good_city].desired_time;
            }
        }
    }

    if(debug)
        cout << "Done with the while loop, will print now" << endl;

    for(int i = 1; i <= n_days; i++){
        for(auto &it: final_ans[i]){
            cout << it << " ";
        }
        cout << endl;
    }


    return 0;
}