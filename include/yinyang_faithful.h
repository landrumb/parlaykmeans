#include <tuple>


template<typename T>
struct yinyang_faithful {

    struct point {

        parlay::slice<T*, T*> coordinates; // the coordinates of the point
        size_t best; // the index of the best center for the point
        size_t id;
        float ub;
        float lb;

    
        point() : best(-1), coordinates(nullptr, nullptr), 
        ub(std::numeric_limits<float>::max()), lb(-1) {
        }

        point(size_t chosen, parlay::slice<T*,T*> coordinates) : best(chosen), 
        coordinates(coordinates.begin(),coordinates.end()), 
        ub(std::numeric_limits<float>::max()), lb(-1) {
        
        }

        point(size_t chosen, parlay::slice<T*,T*> coordinates, float ub, 
        float lb) : best(chosen), 
        coordinates(coordinates.begin(),coordinates.end()), 
        ub(ub), lb(lb) {
        
        }



    };

    struct center {
        size_t id; // a unique (hopefully) identifier for the center
        parlay::sequence<float> coordinates; // the pointer to coordinates of the center
    
        center(size_t id, parlay::sequence<float> coordinates) : id(id) {
        
            this->coordinates = coordinates;
        }

        center() : id(-1) {
        
        }

    };



    void cluster(T* v, size_t n, size_t d, size_t k, float* c, size_t* asg, 
    Distance& D, size_t max_iter, double epsilon) {
        return;
    }
};