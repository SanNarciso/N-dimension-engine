//#include <iostream>
//#include <cmath>
//#include <vector>
//#include <algorithm>
//#include <string>
//#include<random>
//#define SCREEN_WIDTH 1920
//#define SCREEN_HEIGHT 1080
//using namespace std;
//
//uniform_real_distribution<double> unif(1,100);
//default_random_engine re;
//
//
//class Matrix{
//public:
//    int n, m;
//    bool transposed = false;
//
//    vector<vector<double>> matrix;
//
//    Matrix(int rows = 3, int cols = 3, vector<vector<double>> vals = {{1,0,0}, {0,1,0}, {0,0,1}}){
//        n = rows;
//        m = cols;
//        matrix = vector<vector<double>> (n, vector<double> (m, 0));
//        for (int i = 0; i < n; ++i){
//            for (int j = 0; j < m; ++j){
//                matrix[i][j] = vals[i][j];
//            }
//        }
//    }
//
//    virtual void print(){
//        for (int i = 0; i < n; ++i){
//            for (int j = 0; j < m; ++j){
//                cout << matrix[i][j] << ' ';
//            }
//            cout << endl;
//        }
//    }
//
//    Matrix operator+ (Matrix matr){
//        if (matr.n != n || matr.m != m) {
//            cout << "operation is not allowed" << endl;
//            return Matrix(n, m);
//        }
//
//        Matrix result_matrix = *(new Matrix(matr.n, matr.m));
//        for (int i = 0; i < n; ++i){
//            for (int j = 0; j < m; ++j){
//                result_matrix.matrix[i][j] = (*this).matrix[i][j] + matr.matrix[i][j];
//            }
//        }
//        return result_matrix;
//    }
//
//    Matrix operator* (double alpha) {
//        Matrix result_matrix = *(new Matrix(n, m));
//        result_matrix.matrix = matrix;
//        for (int i = 0; i < n; ++i){
//            for (int j = 0; j < m; ++j){
//                result_matrix.matrix[i][j] *= alpha;
//            }
//        }
//        return result_matrix;
//    }
//
//
//    Matrix operator* (Matrix mat){
//        if (m != mat.n){
//            cout << "operation in not allowed" << endl;
//            return Matrix(n, m);
//        }
//
//        Matrix result_matrix(n, mat.m);
//        for (int i = 0; i < n; ++i){
//            for (int j = 0; j < mat.m; ++j){
//                double sum = 0.0;
//                for (int k = 0; k < m; ++k){
//                    sum += (*this).matrix[i][k]*mat.matrix[k][j];
//                }
//                result_matrix.matrix[i][j] = sum;
//            }
//        }
//        return result_matrix;
//    }
//
//    Matrix operator- (Matrix mat){
//        if (mat.n != n || mat.m != m){
//            cout << "operation is not allowed" << endl;
//            return Matrix(mat.n, mat.m);
//        }
//        Matrix result_matrix = *this + mat*(-1);
//        return result_matrix;
//    }
//
//    Matrix operator/ (double alpha){
//        if (alpha == 0){
//            cout << "operation is not allowed" << endl;
//            return Matrix(n,m);
//        }
//        Matrix result_matrix = *this * (1/alpha);
//        return result_matrix;
//    }
//
//    bool operator== (Matrix mat){
//        if (mat.n != n || mat.m != m){
//            cout << "operation is not allowed" << endl;
//            return false;
//        }
//        for (int i = 0; i < n; ++i){
//            for(int j = 0; j < m; ++j){
//                if (matrix[i][j] != mat.matrix[i][j]) return false;
//            }
//        }
//        return true;
//    }
//
//    Matrix sub_matrix(Matrix mat, int i_row, int j_col){
//        Matrix result(mat.n-1, mat.m-1);
//        for (int i = 0; i < mat.n; ++i){
//            for (int j = 0; j < mat.m; ++j){
//                if (i < i_row && j < j_col){
//                    result.matrix[i][j] = mat.matrix[i][j];
//                }
//                else if (i < i_row && j > j_col){
//                    result.matrix[i][j-1] = mat.matrix[i][j];
//                }
//                else if (i > i_row && j < j_col){
//                    result.matrix[i-1][j] = mat.matrix[i][j];
//                }
//                else if (i > i_row && j > j_col){
//                    result.matrix[i-1][j-1] = mat.matrix[i][j];
//                }
//            }
//        }
//        return result;
//    }
//
//     double det() {
//         if (n != m) {
//             cout << "operation DET is not allowed" << endl;
//             return 0;
//         }
//         if (n == 1) return (*this).matrix[0][0];
//         double deter = 0;
//         for (int COL = 0; COL < (*this).matrix.size(); ++COL){
//             Matrix sub_mat = sub_matrix(*this, 0, COL);
//             deter += pow(-1, COL) * matrix[0][COL] * sub_mat.det();
//         }
//         return deter;
//     }
//
//     bool exists_inverse(){
//         if ((*this).det() == 0 || n != m) {
//             cout << "inverse does not exist" << endl;
//             return false;
//         }
//         return true;
//    }
//
//     Matrix inverse(){
//         if (!(*this).exists_inverse()) return (*this);
//
//        Matrix copy = *this;
//        Matrix identity(n, n);
//
//        // upper triangle
//         for (int i = 0; i < n; i++) {
//             if (copy.matrix[i][i] == 0) {
//                 for (int j = i + 1; j < n; j++) {
//                     if (copy.matrix[j][i] != 0) {
//                         swap(copy.matrix[i], copy.matrix[j]);
//                         swap(identity.matrix[i], identity.matrix[j]);
//                         break;
//                     }
//                 }
//             }
//             double factor = copy.matrix[i][i];
//             for (int j = 0; j < n; j++) {
//                 copy.matrix[i][j] /= factor;
//                 identity.matrix[i][j] /= factor;
//             }
//             for (int j = i + 1; j < n; j++) {
//                 double factor = copy.matrix[j][i];
//                 for (int k = 0; k < n; k++) {
//                     copy.matrix[j][k] -= factor * copy.matrix[i][k];
//                     identity.matrix[j][k] -= factor * identity.matrix[i][k];
//                 }
//             }
//         }
//
//         // diagonal matrix
//         for (int i = n - 1; i >= 0; i--) {
//             for (int j = i - 1; j >= 0; j--) {
//                 double factor = copy.matrix[j][i];
//                 for (int k = 0; k < n; k++) {
//                     copy.matrix[j][k] -= factor * copy.matrix[i][k];
//                     identity.matrix[j][k] -= factor * identity.matrix[i][k];
//                 }
//             }
//         }
//         return identity;
//    }
//
//    Matrix operator/ (Matrix M){
//        if (!(M.exists_inverse())) return (*this);
//        return (*this)*M.inverse();
//    }
//
//    Matrix transpose(){
//        Matrix res(m, n);
//        for (int i = 0; i < m; ++i){
//            for (int j = 0; j < n; ++j){
//                res.matrix[i][j] = matrix[j][i];
//            }
//        }
//        !(*this).transposed ? (*this).transposed = true : (*this).transposed = false;
//        return res;
//    }
//
//
//    // HERE YOU CAN CHANGE THE NORM OF MATRIX
//    double norm(){
//        double res = 0;
//        for (int i = 0; i < n; ++i){
//            for (int j = 0; j < m; ++j){
//                res += pow(matrix[i][j], 2);
//            }
//        }
//        return sqrt(res);
//    }
//
//    static Matrix rotation_x(double alpha){
//        return Matrix(3,3, {{1,0,0},{0, cos(alpha), -sin(alpha)}, {0, sin(alpha), cos(alpha)}});
//    }
//    static Matrix rotation_y(double alpha){
//        return Matrix(3,3, {{cos(alpha), 0, sin(alpha)},{0, 1, 0}, {-sin(alpha), 0, cos(alpha)}});
//    }
//    static Matrix rotation_z(double alpha){
//        return Matrix(3,3, {{cos(alpha),-sin(alpha),0},{sin(alpha), cos(alpha), 0}, {0, 0, 1}});
//    }
//};
//
//class Vector : public Matrix{
//
//public:
//    int dim;
//    double length;
//    vector<double> coordinates;
//    vector<double> coordinates_in_basis;
//
//    /* NOTE: coordinates and coordinates_in_basis are not necessary the same size:
//     *      vector with coordinates (a_1, a_2, a_3) in basis {(1,0,0, 0), (0,1,0,0), (0,0,1,0), (0,0,0,1)} has coordinates_in_basis = (a_1, a_2, a_3, 0)
//     */
//
//    Vector(int size = 3) : Matrix(size, 1){
//        dim = size;
//        for (int i = 0; i < size; ++i){
//            coordinates.push_back(matrix[i][0]);
//        }
//    }
//
//    Vector (vector<double> coords) : Matrix(coords.size(), 1){
//        dim = coords.size();
//        for (int i = 0; i < coords.size(); ++i){
//            coordinates.push_back(coords[i]);
//            matrix[i][0] = coords[i];
//        }
//    }
//
//    void change_coordinates(vector<double> coord, vector<int> indicies){
//        if (coord.size() != indicies.size()) throw invalid_argument("Invalid transformation of coordinates");
//        for (int i = 0; i < coord.size(); ++i) {
//            coordinates[indicies[i]] = coord[i];
//            matrix[indicies[i]][0] = coord[i];
//        }
//
//    }
//
//    Vector operator+ (Vector v){
//        Matrix result = Matrix::operator+(v);
//        return get_column(result, 0);
//    }
//
//    Vector operator* (double a){
//        Matrix result = Matrix::operator*(a);
//        return get_column(result, 0);
//    }
//
//    Vector operator- (Vector v){
//        Matrix result = Matrix::operator-(v);
//        return get_column(result, 0);
//    }
//
//    Vector operator/ (double a){
//        Matrix result = Matrix::operator/(a);
//        return get_column(result, 0);
//    }
//
//
//    Matrix operator*(Matrix M){
//        Matrix result = Matrix::operator*(M);
//        return result;
//    }
//
//    static Vector forAll(int sz){
//        Vector v(sz);
//        for (int i = 0; i < sz; ++i){
//            v.coordinates[i] = unif(re);
//            v.matrix[i][0] = v.coordinates[i];
//        }
//        return v;
//    }
//
//    bool isZero(){
//        for (auto i: coordinates){
//            if (i != 0) return false;
//        }
//        return true;
//    }
//
//    Vector operator^ (Vector v){
//        if ((*this).coordinates.size() != v.coordinates.size() || v.coordinates.size() != 3) {
//            throw invalid_argument("Input vectors must have the same size and be of length 3");
//        }
//        Vector res(3);
//        res.coordinates[0] = (*this).coordinates[1] * v.coordinates[2] - (*this).coordinates[2] * v.coordinates[1];
//        res.coordinates[0] = (*this).coordinates[2] * v.coordinates[0] - (*this).coordinates[0] * v.coordinates[2];
//        res.coordinates[0] = (*this).coordinates[0] * v.coordinates[1] - (*this).coordinates[1] * v.coordinates[0];
//        return res;
//    }
//private:
//    Vector get_column(Matrix M, int col){
//        Vector result(M.matrix.size());
//        for (int i = 0; i < M.matrix.size(); ++i){
//            result.coordinates[i] = M.matrix[i][col];
//            result.matrix[i][col] = M.matrix[i][col];
//        }
//        return result;
//    }
//    Vector get_row(Matrix M, int row){
//        Vector result(M.matrix.size());
//        for (int i = 0; i < M.matrix.size(); ++i){
//            result.coordinates[i] = M.matrix[row][i];
//            result.matrix[row][i] = M.matrix[row][i];
//        }
//        return result;
//    }
//};
//
//
//class BilinearForm{
//public:
////    double result = 0;
//    bool symmetric = false;
//    bool positive_defined = false;
//    Matrix B_mat;
//    Matrix B_mat_in_basis;                      // ONLY WHEN BASIS APPEARS
//
//    double evaluate(Vector v1, Vector v2, Matrix M){
//        Matrix res = v1.transpose() * M * v2;
//        return res.matrix[0][0];
//    }
//
//    BilinearForm(Vector v1 = Vector(), Vector v2 = Vector(), Matrix M = Matrix(), int dim = 3){
//        if (M.n != M.m || dim != M.n){
//            cout << "incorrect matrix or dim for building bilinear from" << endl;
//            return;
//        }
//
//        //exceptions at the beginning of the method
//
//        Vector v3 = Vector::forAll(dim);
//        double alpha = unif(re);
//        double beta = unif(re);
//
//        double temp0 = ceil(evaluate(v1+v3,v2, M));
//        double temp1 = ceil(evaluate(v1,v2, M) + evaluate(v3,v2, M));
//        double temp2 = ceil(evaluate(v1,v2 + v3, M));
//        double temp3 = ceil(evaluate(v1,v2, M) + evaluate(v1,v3, M));
//        double temp4 = ceil(evaluate(v1*alpha,v2*beta, M));
//        double temp5 = ceil(alpha*beta*evaluate(v1,v2, M));
//        if (temp0 == temp1 && temp2 == temp3 && temp4 == temp5) {
//            //result = evaluate(v1, v2, M);
//            B_mat = M;
//        }
//        else{
//            cout << "Bad matrix for bilinear form: does not satisfy axioms" << endl;
//        }
//
//        if (evaluate(v1,v2,M) == evaluate(v2,v1,M)) symmetric = true;
//
//        if ((*this).symmetric && evaluate(v3, v3, M) > 0){
//            positive_defined = true;
//        } else cout << "Matrix does not form scalar product" << endl;
//    }
//};
//
//
//
//class EuclidSpace{
//public:
//    vector<Vector> basis_vectors;
//    Matrix basis_vectors_matrix;
//    BilinearForm scalar_prod;
//    int dimension_of_vector_space = 0;
//    Matrix Gramm;
//
//    EuclidSpace(vector<Vector> elements, BilinearForm sp){
//        dimension_of_vector_space = elements.size();
//        basis_vectors_matrix = Matrix(dimension_of_vector_space, dimension_of_vector_space);
//        Matrix temp(elements.size(), elements.size());
//        for (int i = 0; i < elements.size(); ++i){
//            for (int j = 0; j < elements.size(); ++j){
//                temp.matrix[i][j] = elements[i].coordinates[j];
//            }
//        }
//        if (!temp.exists_inverse()){
//            cout << "Vectors are linearly dependent" << endl;
//            return;
//        }
//
//        for (int i = 0; i < elements.size(); ++i){
//            if (elements[i].dim != dimension_of_vector_space){
//                cout << "Elements do not form a basis: bad dimensions" << endl;
//                return;
//            }
//            else{
//                elements[i].coordinates_in_basis = vector<double> (dimension_of_vector_space, 0);
//                elements[i].coordinates_in_basis[i] = 1;
//                basis_vectors.push_back(elements[i]);
//            }
//        }
//
//        for (int i = 0; i < dimension_of_vector_space; ++i){
//            for (int j = 0; j < dimension_of_vector_space; ++j){
//                basis_vectors_matrix.matrix[i][j] = basis_vectors[j].coordinates[i];
//            }
//        }
//
//        if (sp.symmetric && sp.positive_defined){
//            scalar_prod = sp;
//            for (int i = 0; i < elements.size(); ++i){
//                for (int j = 0; j < elements.size(); ++j){
//                    scalar_prod.B_mat_in_basis.matrix[i][j] = scalar_prod.evaluate(elements[i], elements[j], sp.B_mat);
//                }
//            }
//            Gramm = scalar_prod.B_mat_in_basis;
//        } else cout << "Current bilinear form does not form a scalar product";
//    }
//
//    void get_basis_coordinates(Vector &v){
//        if (v.dim > dimension_of_vector_space){
//            throw invalid_argument("Impossible to decompose vector in basis: dim v > dim E where E is current euclid space");
//        }
//
//        if (v.dim == dimension_of_vector_space) v.coordinates_in_basis = (((basis_vectors_matrix.inverse()) * v).transpose()).matrix[0];
//        else{
//            Vector copy(dimension_of_vector_space);
//            for (int i = 0; i < v.dim; ++i){
//                copy.coordinates[i] = v.coordinates[i];
//                copy.matrix[i][0] = v.matrix[i][0];
//            }
//            v.coordinates_in_basis = (((basis_vectors_matrix.inverse()) * copy).transpose()).matrix[0];
//        }
//
//        Vector tmp = Vector(v.coordinates_in_basis);
//        v.length = sqrt((tmp.transpose() * Gramm * tmp).matrix[0][0]);
//    }
//};
//
//class Point{
//public:
//    vector<double> coordinates;
//    Point(Vector v){
//        coordinates = v.coordinates_in_basis;
//    }
//
//    Point operator+(Vector v) {
//        Point res(Vector(coordinates.size()));             // dim v == dim v_in_basis ? In general, no.
//        for (int i = 0; i < coordinates.size(); ++i) {
//            res.coordinates[i] += v.coordinates_in_basis[i];
//        }
//        return res;
//    }
//};
//
//
//class CoordinateSystem{
//    /*
//     * What changes when coordinate system appears?
//     * Basis is enough, isn't it?
//     */
//public:
//};
///*
//int main(){
//
////    Matrix t = Matrix(3,3);
////    Matrix t1 = Matrix(3, 3);
//
//    //t.matrix[0][0] = 10;
//    //cout << t.det() << endl;
//    //t.print();
//    //(t.inverse()).print();
//    //(t*t.inverse()).print();
//
//    Vector v(3);
//    //v.change_coordinates({5,5,5}, {0,1,2});
//
//    //cout << "v: " << endl;
//    //v.print();
//    //cout << endl;
//
////    Vector v1(3);
////    v1.change_coordinates({1,1,1}, {0,1,2});
//
//    //v1.print(); cout << endl;
//    //Vector sum = v+v1;
//    //sum.print();
//    //(v/3).print();
//    //cout << typeid(v+v1).name() << endl;
//    //Vector q = Vector::forAll(3);
//    //Vector q1 = Vector::forAll(3);
//    //q.print();
//    //q1.print();
//
////    BilinearForm B(v,v1,t, t.n);
////    //cout << B.result << ' ' << B.symmetric << endl;
////    EuclidSpace E({Vector({1,0,0}), Vector({0,2,0}), Vector({0,0,1})}, B);
////    E.scalar_prod.B_mat.print();
////    E.Gramm.print();
////    Vector test({1,1,0});
////    E.get_basis_coordinates(test);
////    test.print();
////    for (auto i : test.coordinates_in_basis) cout << i << endl;
////    cout << test.length;
//
//    //Vector v3 = Vector::forAll(3);
//    //v.print();
//    //v1.print();
//    //Vector s = v+v1;
//    //s.print();
//    //v3.print();
//    //((v.transpose()) * t * v1).print();
//    //cout << evaluate(v3,v,t);
//    //v.print();
////    Matrix m2 = v.transpose() * v;  // m2.type == Matrix because when
////    Matrix m3 = v * v.transpose();  // this MUST be a Matrix instance because ( n x 1) x (1 x n) -> n x n
////    m2.print();
//    Matrix::rotation_x(M_PI/2).print();
//    (Matrix::rotation_y(M_PI/3) * v).print();
//    cout << cos(3.14/2);
//}*/
//
