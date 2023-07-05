/* 
    A struct which wraps a flat array representing points to facilitate access to the points' coordinates.
 */

template <typename T>
struct Points {
    T* points;
    size_t num_points;
    size_t dim;

    size_t width = sizeof(T); // width of each dimension in bytes

    Points(T* points, size_t num_points, size_t dim) : points(points), num_points(num_points), dim(dim) {}

    T* operator[](size_t i) {
        return points + (i * dim) * width;
    }
};

// TODO add points struct that handles mmapped fvec or bvec files