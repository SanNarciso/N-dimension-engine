#include <gtest/gtest.h>
#include "../../lib/Engine/RayCast.h"
class MatrixTest : public ::testing::Test {
public:

    Matrix matrix_1;
    Matrix matrix_2;
    Matrix matrix_3;
};


TEST_F(MatrixTest, TestMatrixAddition) {
    matrix_1 = Matrix(3, 3, {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    matrix_2 = Matrix(3, 3, {{9, 8, 7}, {6, 5, 4}, {3, 2, 1}});
    matrix_3 = matrix_1 + matrix_2;
    Matrix mat = Matrix(3, 3, {{10,10,10}, {10,10,10}, {10,10,10}});
    ASSERT_EQ(matrix_3, mat);
}

TEST_F(MatrixTest, TestMatrixOutOfBounds) {
    matrix_1 = Matrix(3, 3, {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    matrix_2 = Matrix(2, 2, {{9, 8}, {7, 6}});
//    matrix_3 = matrix_1 + matrix_2;
    EXPECT_THROW(matrix_1 + matrix_2, MatrixException);
}

TEST(Matrix, MultiplicationWithScalar) {
    Matrix mat(2, 2, {{1, 2}, {3, 4}});
    Matrix result = mat * 2.0;

    EXPECT_EQ(result.matrix[0][0], 2.0);
    EXPECT_EQ(result.matrix[0][1], 4.0);
    EXPECT_EQ(result.matrix[1][0], 6.0);
    EXPECT_EQ(result.matrix[1][1], 8.0);
}

TEST(Matrix, MultiplicationWithAnotherMatrix) {
    Matrix mat1(2, 3, {{1, 2, 3}, {4, 5, 6}});
    Matrix mat2(3, 2, {{7, 8}, {9, 10}, {11, 12}});
    Matrix result = mat1 * mat2;

    EXPECT_EQ(result.matrix[0][0], 58.0);
    EXPECT_EQ(result.matrix[0][1], 64.0);
    EXPECT_EQ(result.matrix[1][0], 139.0);
    EXPECT_EQ(result.matrix[1][1], 154.0);
}

TEST(Matrix, InvalidMultiplicationWithAnotherMatrix) {
    Matrix mat1(2, 3, {{1, 2, 3}, {4, 5, 6}});
    Matrix mat2(2, 2, {{7, 8}, {9, 10}});

    EXPECT_THROW(mat1*mat2, MatrixException);
}

TEST(MatrixTests, Subtraction) {
    Matrix m1(2, 2, {{1, 2}, {3, 4}});
    Matrix m2(2, 2, {{5, 6}, {7, 8}});
    Matrix expected(2, 2, {{-4, -4}, {-4, -4}});

    Matrix result = m1 - m2;

    EXPECT_EQ(result, expected);
}

TEST(MatrixTests, Division) {
    Matrix m1(2, 2, {{1, 2}, {3, 4}});
    double alpha = 2.0;
    Matrix expected(2, 2, {{0.5, 1}, {1.5, 2}});

    Matrix result = m1 / alpha;

    EXPECT_EQ(result, expected);
}

TEST(MatrixTests, DivisionByZero){
    EXPECT_THROW(Matrix(2, 2, {{1, 2}, {3, 4}}) / 0, MatrixException);
}

TEST(MatrixTests, Equality) {
    Matrix m1(2, 2, {{1, 2}, {3, 4}});
    Matrix m2(2, 2, {{1, 2}, {3, 4}});
    Matrix m3(2, 2, {{1, 2}, {4, 3}});

    EXPECT_EQ(m1, m2);
    EXPECT_NE(m1, m3);
}

TEST(SubMatrix, SubMatrixReturnsCorrectSubmatrix) {
    Matrix mat(3, 3);
    mat.matrix = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    Matrix submat = mat.sub_matrix(1, 1);
    EXPECT_EQ(submat.matrix, (vector<vector<double>>{{1, 3}, {7, 9}}));
}

TEST(SubMatrix, SubMatrixReturnsIdentityMatrixForInvalidIndices) {
    Matrix mat(3, 3);
    mat.matrix = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    EXPECT_THROW(mat.sub_matrix(4,4), MatrixException);
}

TEST(Det, DetReturnsCorrectDeterminant) {
    Matrix mat(3, 3);
    mat.matrix = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    double deter = mat.det();
    EXPECT_EQ(deter, 0);
}

TEST(Det, DetReturnsCorrectDeterminantFor2x2Matrix) {
    Matrix mat(2, 2);
    mat.matrix = {{1, 2}, {3, 4}};
    double deter = mat.det();
    EXPECT_EQ(deter, -2);
}

TEST(Det, DetReturnsCorrectDeterminantFor3x3Matrix) {
    Matrix mat(3, 3);
    mat.matrix = {{2, 3, 4}, {5, 6, 7}, {8, 9, 1}};
    double deter = mat.det();
    EXPECT_EQ(deter, 27);
}



TEST_F(MatrixTest, OperatorDivide) {
    Matrix mat1(2, 2);
    mat1.matrix = {{2, 4}, {6, 8}};
    Matrix mat2(2, 2);
    mat2.matrix = {{1, 2}, {3, 4}};
    Matrix result = mat1 / mat2;
    Matrix expected_result(2, 2);
    expected_result.matrix = {{2, 0}, {0, 2}};
    ASSERT_TRUE(result == expected_result);

    // Test for a non-invertible matrix
    Matrix mat3(2, 2);
    mat3.matrix = {{1, 2}, {2, 4}};
    EXPECT_THROW(mat1/mat3, MatrixException);

}

TEST_F(MatrixTest, Transpose) {
    Matrix mat(3, 2);
    mat.matrix = {{1, 2}, {3, 4}, {5, 6}};
    Matrix result = mat.transpose();
    Matrix expected_result(2, 3);
    expected_result.matrix = {{1, 3, 5}, {2, 4, 6}};
    ASSERT_TRUE(result == expected_result);
    ASSERT_TRUE(mat.transposed);
    result = result.transpose();
    ASSERT_TRUE(result == mat);
}

TEST(VectorTest, Rotation2d){
    Vector v(2);
    Matrix rotated = Matrix::Plane_rotate(M_PI/3)*v;
    rotated.print();
    EXPECT_EQ(rotated.matrix[0][0], 0.5);
//    EXPECT_EQ(rotated.matrix[1][0], 0.866025);
}

TEST(VectorTest, Rotation_x_y_z){

}


TEST(VectorTest, DefaultConstructor){
    Vector v;
    EXPECT_EQ(v.dim, 3);
    EXPECT_EQ(v.coordinates.size(), 3);

    EXPECT_EQ(v.coordinates[0], 1);
    EXPECT_EQ(v.matrix[0][0], 1);
    EXPECT_EQ(v.coordinates[1], 0.0);
    EXPECT_EQ(v.matrix[1][0], 0.0);
    EXPECT_EQ(v.coordinates[2], 0.0);
    EXPECT_EQ(v.matrix[2][0], 0.0);
}

TEST(VectorTest, ParameterizedConstructor){
    vector<double> coords = {1.0, 2.0, 3.0};
    Vector v(coords);
    EXPECT_EQ(v.dim, 3);
    EXPECT_EQ(v.coordinates.size(), 3);
    for(int i = 0; i < 3; i++){
        EXPECT_EQ(v.coordinates[i], coords[i]);
        EXPECT_EQ(v.matrix[i][0], coords[i]);
    }
}


TEST(VectorTest, CrossProduct){
    vector<double> coords1 = {1.0, 2.0, 3.0};
    Vector v1(coords1);
    vector<double> coords2 = {4.0, 5.0, 6.0};
    Vector v2(coords2);
    vector<double> expected_coords = {-3.0, 6.0, -3.0};
    Vector expected(expected_coords);
    Vector result = v1^v2;
    EXPECT_EQ(result.dim, 3);
    EXPECT_EQ(result.coordinates.size(), 3);
    for(int i = 0; i < 3; i++){
        EXPECT_DOUBLE_EQ(result.coordinates[i], expected.coordinates[i]);
        EXPECT_DOUBLE_EQ(result.matrix[i][0], expected.matrix[i][0]);
    }
}


//TEST(EntityTest, EntitySetPropertyTest){
//    Entities::Entity entity;
//    entity.set_property(map<string, int> {{"entity_property_function_test", 0}});
//    map<string, int>::iterator it = entity.properties.begin();
//    EXPECT_EQ(it->first, "entity_property_function_test");
//    EXPECT_EQ(it->second, 0);
//}

//TEST(EntityTest, EntityGetPropertyTest){
//    Entities::Entity entity;
//    entity.set_property(map<string, int> {{"entity_property_function_test", 0}});
//    map<string, int>::iterator it = entity.properties.begin();
//    EXPECT_EQ(entity.get_property("entity_property_function_test"), 0);
//}

//TEST(EntityTest, EntityRemovePropertyTest){
//    Entities::Entity entity;
//    entity.set_property(map<string, int> {{"entity_property_function_test", 0}});
//    map<string, int>::iterator it = entity.properties.begin();
//    entity.remove_property("entity_property_function_test");
//    EXPECT_EQ(entity.properties.empty(), true);
//}


TEST(RAYCASTING, RayCastConstruct){
    /*
     * HERE IS A PROBLEM OCCURED:
     * each time we create a coordinate system, we need to get its basis coordinates to assign a value to field length in this Vector class instance
     */
    RayCast ray;
    cout << ray.cam.fov;
}


//int main(){
//    testing::InitGoogleTest();
//    RUN_ALL_TESTS();
//}
