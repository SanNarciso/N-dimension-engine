#pragma once
#include "lib/Engine/Entities.h"

int main(){
        /*
     * To create a scene, we need to:
     * 1) create an Ellipsoid (for example), here it is added to EntitiesList
     * 2) create a Camera
     * 3) create a Console
     * 4) draw
     */

    Entities::Game::HyperEllipsoid HE(Point(Vector({1000,1000,1000})), Vector({0.1,0.1,0.1}), {Vector({0.3, 0, 0}), Vector({0, 0.3, 0}), Vector({0, 0, 0.3})});
    //Entities::Entity* ent = &HE;

    Entities::Game::Camera camera(60, INT_MAX);
    auto q = Entities::Entities_List;
    camera.position = Point(Vector({0,0,0}));
    camera.direction = Vector({1,1,1});
    camera.cs.basis.get_basis_coordinates(camera.direction);
    vector<vector<Entities::Ray>> Q = camera.get_rays_matrix(40, 40);
    Entities::Game::Canvas console(camera, 40,40);
    console.update(HE);
    //cout << HP.intersection_distance(Entities::Ray(HP.cs, Point(Vector({0,0,0})), Vector({1,2,1})));
    console.distances.print(); // trouble's here
    //console.drawObjects(HE);
}