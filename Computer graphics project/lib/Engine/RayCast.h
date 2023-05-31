#ifndef ENGINE_RAYCAST_H
#define ENGINE_RAYCAST_H

#endif //ENGINE_RAYCAST_H
#include "Entities.h"

class RayCast{
public:
    int n,m;
    vector<vector<Vector>> A;
    Entities::Game::Camera cam;
    Vector directive;
    double alpha, beta;
    RayCast(int nn = 16, int mm = 9){
        n = nn;
        m = mm;
        A = vector<vector<Vector>> (n, vector<Vector> (m));
        alpha = cam.fov;
        beta = cam.vfov;
        directive = cam.direction;
        double d_alpha = alpha/n;
        double d_beta = beta/m;
        vector<double> alpha_i(n, 0); vector<double> beta_j(m,0);
        for (int i = 0; i < n; ++i){
            alpha_i[i] = d_alpha*i - alpha/2;
        }
        for (int j = 0; j < m; ++j){
            beta_j[j] = d_beta*j - beta/2;
        }

        for (int i = 0; i < n; ++i){
            for (int j = 0; j < m; ++j){
                Vector v_i_j((((Matrix::rotation_z(alpha_i[i])*Matrix::rotation_y(beta_j[j]))*directive).transpose()).matrix[0]);
                Vector _v_i_j = v_i_j*(pow(directive.length, 2)/((directive.transpose()*cam.cs.basis.Gramm*v_i_j).matrix[0][0]));
                A[i][j] = _v_i_j;
            }
        }
    }
};