#include <iostream>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <array>
#include <cmath>
#include <random>

using namespace std;

/**
 * Class for representing a point. coordinate_type must be a numeric type.
 */
template<typename coordinate_type, size_t dimensions>
class point {
public:
    point(array<coordinate_type, dimensions> c) : coords_(c) {}
    point(initializer_list<coordinate_type> list) {
        size_t n = min(dimensions, list.size());
        copy_n(list.begin(), n, coords_.begin());
    }
    /**
     * Returns the coordinate in the given dimension.
     *
     * @param index dimension index (zero based)
     * @return coordinate in the given dimension
     */
    coordinate_type get(size_t index) const {
        return coords_[index];
    }
    /**
     * Returns the distance squared from this point to another
     * point.
     *
     * @param pt another point
     * @return distance squared from this point to the other point
     */
    double distance(const point& pt) const {
        double dist = 0;
        for (size_t i = 0; i < dimensions; ++i) {
            double d = get(i) - pt.get(i);
            dist += d * d;
        }
        return dist;
    }
private:
    array<coordinate_type, dimensions> coords_;
};

template<typename coordinate_type, size_t dimensions>
ostream& operator<<(ostream& out, const point<coordinate_type, dimensions>& pt) {
    out << '(';
    for (size_t i = 0; i < dimensions; ++i) {
        if (i > 0)
            out << ", ";
        out << pt.get(i);
    }
    out << ')';
    return out;
}

/**
 * C++ k-d tree implementation, based on the C version at rosettacode.org.
 */
template<typename coordinate_type, size_t dimensions>
class kdtree {
public:
    typedef point<coordinate_type, dimensions> point_type;

private:
    struct node {
        node(const point_type& pt) : point_(pt), left_(nullptr), right_(nullptr) {}
        coordinate_type get(size_t index) const {
            return point_.get(index);
        }
        double distance(const point_type& pt) const {
            return point_.distance(pt);
        }
        point_type point_;
        node* left_;
        node* right_;
    };
    node* root_ = nullptr;
    node* best_ = nullptr;
    double best_dist_ = 0;
    size_t visited_ = 0;
    vector<node> nodes_;

    struct node_cmp {
        node_cmp(size_t index) : index_(index) {}
        bool operator()(const node& n1, const node& n2) const {
            return n1.point_.get(index_) < n2.point_.get(index_);
        }
        size_t index_;
    };

    node* make_tree(size_t begin, size_t end, size_t index) {
        if (end <= begin)
            return nullptr;
        size_t n = begin + (end - begin)/2;
        auto i = nodes_.begin();
        nth_element(i + begin, i + n, i + end, node_cmp(index));
        index = (index + 1) % dimensions;
        nodes_[n].left_ = make_tree(begin, n, index);
        nodes_[n].right_ = make_tree(n + 1, end, index);
        return &nodes_[n];
    }

    void nearest(node* root, const point_type& point, size_t index) {
        if (root == nullptr)
            return;
        ++visited_;
        double d = root->distance(point);
        if (best_ == nullptr || d < best_dist_) {
            best_dist_ = d;
            best_ = root;
        }
        if (best_dist_ == 0)
            return;
        double dx = root->get(index) - point.get(index);
        index = (index + 1) % dimensions;
        nearest(dx > 0 ? root->left_ : root->right_, point, index);
        if (dx * dx >= best_dist_)
            return;
        nearest(dx > 0 ? root->right_ : root->left_, point, index);
    }

public:
    kdtree(const kdtree&) = delete;
    kdtree& operator=(const kdtree&) = delete;
    /**
     * Constructor taking a pair of iterators. Adds each
     * point in the range [begin, end) to the tree.
     *
     * @param begin start of range
     * @param end end of range
     */
    template<typename iterator>
    kdtree(iterator begin, iterator end) : nodes_(begin, end) {
        root_ = make_tree(0, nodes_.size(), 0);
    }
    
    /**
     * Constructor taking a function object that generates
     * points. The function object will be called n times
     * to populate the tree.
     *
     * @param f function that returns a point
     * @param n number of points to add
     */
    template<typename func>
    kdtree(func&& f, size_t n) {
        nodes_.reserve(n);
        for (size_t i = 0; i < n; ++i)
            nodes_.push_back(f());
        root_ = make_tree(0, nodes_.size(), 0);
    }

    /**
     * Returns true if the tree is empty, false otherwise.
     */
    bool empty() const { return nodes_.empty(); }

    /**
     * Returns the number of nodes visited by the last call
     * to nearest().
     */
    size_t visited() const { return visited_; }

    /**
     * Returns the distance between the input point and return value
     * from the last call to nearest().
     */
    double distance() const { return sqrt(best_dist_); }

    /**
     * Finds the nearest point in the tree to the given point.
     * It is not valid to call this function if the tree is empty.
     *
     * @param pt a point
     * @return the nearest point in the tree to the given point
     */
    const point_type& nearest(const point_type& pt) {
        if (root_ == nullptr)
            throw logic_error("tree is empty");
        best_ = nullptr;
        visited_ = 0;
        best_dist_ = 0;
        nearest(root_, pt, 0);
        return best_->point_;
    }
};

class City{
    public:
        int street, avenue, desired_time, start_hour, end_hour, mid_hour, value, cluster_number;
};

unordered_map<int, vector<City>> make_clusters(vector<City> cities, int n_clusters = 5){
    // Assume n_clusters is fixed for now
    // (TODO: Fix this)
    unordered_map<int, pair<pair<pair<int, int>, int>, vector<City>>> clusters;
    // Hash map: custers
    // Key: cluster_number
    // Value: Pair of cluster center and list of cities
    // Cluster center is represented as 3 integers, each in dimension x, y and z
    // as a pair of pair of ints for x and y and an int for z

    // Assign each city a cluster number from 1 to n_clusters
    for (auto &city: cities){
        city.cluster_number = 0; // Use random library to randomly assign a number from 1 to n_clusters
        clusters[0].second.push_back(city);
    }

    int cost, epsilon, n, c_x, c_y, c_z;
    bool converged = false;
    while(converged){
        for(auto &it: clusters){
            c_x = 0; c_y = 0; c_z = 0;
            n = 0;
            for(auto &city: it.second.second){
                c_x += city.avenue; c_y += city.street; c_z += city.mid_hour;
                n++;
            }

            c_x /= n; c_y /= n; c_z /= n;
        }

        for(auto &city: clusters){

        }
        
    }
    
}

int main(){
    unordered_map<int, vector<City>> days_to_city; // Assume this hash map has a list of cities for each day alredy filled out
    unordered_map<int, unordered_map<int, vector<City>>> days_to_cluster;

    // Time to run k means for each day
    for (auto &it: days_to_city){
        auto &clusters = make_clusters(it.second);
        days_to_cluster[it.first] = clusters;
    }

    // Now for every feasible combination of clusters, select the top few ones
    // TODO

    return 0;
}