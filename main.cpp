#include <iostream>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

class City{
    public:
        int street, avenue, mid_hour;
        int id, desired_time, begin_hour, end_hour, value, cluster_number;

        City() {}
        City(int s, int a, int st, int d, int v) : id(s), avenue(a), street(st), desired_time(d), value(v) {}
};

double get_dist(const pair<pair<int, int>, int> &p1, City &city){
    return (p1.first.first - city.avenue) * (p1.first.first - city.avenue) + (p1.first.second - city.street) * (p1.first.second - city.street) + (p1.second - city.mid_hour) * (p1.second - city.mid_hour);
}

int get_nearest_in(vector<pair<pair<int, int>, int>> &points, City &city){
    // TODO: Improve this!!
    if(empty(points)) return -1;

    double nearest_dist = get_dist(points[0], city);
    int nearest = 0;
    for(int i = 0; i < points.size(); i++){
        if(get_dist(points[i], city) < nearest_dist){
            nearest_dist = get_dist(points[i], city);
            nearest = i;
        }
    }

    return nearest;
}

const unordered_map<int, vector<City>> make_clusters(vector<City> &cities, const int n_clusters = 5){
    // Assume n_clusters is fixed for now
    // (TODO: Fix this)
    unordered_map<int, vector<City>> clusters;
    // Hash map: custers
    // Key: cluster_number
    // Value: Pair of cluster center and list of cities
    // Cluster center is represented as 3 integers, each in dimension x, y and z
    // as a pair of pair of ints for x and y and an int for z

    // Assign each city a cluster number from 1 to n_clusters
    for (auto &city: cities){
        city.cluster_number = 0; // TODO: Use random library to randomly assign a number from 1 to n_clusters
        clusters[0].push_back(city);
    }

    vector<pair<pair<int, int>, int>> points;
    points.reserve(n_clusters);

    int cost, epsilon, n;
    double c_x, c_y, c_z;
    bool converged = false;
    while(!converged){
        converged = true;
        for(auto &it: clusters){
            c_x = 0; c_y = 0; c_z = 0;
            n = 0;
            for(auto &city: it.second){
                c_x += city.avenue; c_y += city.street; c_z += city.mid_hour;
                n++;
            }

            c_x /= n; c_y /= n; c_z /= n;
            points.push_back(make_pair(make_pair(c_x, c_y), c_z));
        }

        for(auto &cluster: clusters){
            for(auto &city: cluster.second){
                int new_cluster = get_nearest_in(points, city);
                
                if(new_cluster != city.cluster_number){
                    converged = false;
                }

                city.cluster_number = new_cluster;
            }   
        }
    }

    // TODO: Calculate entropy here, and the choose n_clusters such that entropy is minimum

    return clusters;
}

int main(){
    string line;
    ifstream infile("input.txt");

    unordered_map<int, City> cities;
    unordered_map<int, vector<City>> days_to_city;
    
    getline(infile, line);
    while (getline(infile, line)){
        istringstream iss(line);
        int s, a, st, d, v;
        if (!(iss >> s >> a >> st >> d >> v)) { break; }
        cities[s] = City(s, a, st, d, v);
    }

    getline(infile, line);
    while (getline(infile, line)){
        istringstream iss(line);
        int s, d, b, e;
        if (!(iss >> s >> d >> b >> e)) { break; }
        
        City c = City(cities[s]);
        
        c.begin_hour = b;
        c.end_hour = e;
        c.mid_hour = (b + e) / 2;

        days_to_city[d].push_back(c);
    }

    unordered_map<int, unordered_map<int, vector<City>>> days_to_cluster;

    // Time to run k means for each day
    for (auto &it: days_to_city){
        days_to_cluster[it.first] = make_clusters(it.second);
    }

    // Now for every feasible combination of clusters, select the top few ones
    // TODO

    return 0;
}