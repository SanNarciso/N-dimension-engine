#include "lib/Engine/Entities.h"

int main(){
        /*
     * To create a scene, we need to:
     * 1) create an Ellipsoid (for example), here it is added to EntitiesList
     * 2) create a Camera
     * 3) create a Console
     * 4) draw
     */

    Entities::Game::HyperEllipsoid HE(Point(Vector({2,2,0})), Vector({1,1,1}), {Vector({0.5, 0, 0}), Vector({0, 0.5, 0}), Vector({0, 0, 0.5})});
    //Entities::Entity* ent = &HE;s
    Entities::Game::Camera camera(60, 20);
    auto q = Entities::Entities_List;
    camera.position = Point(Vector({-0.05,-0.05,0}));
    camera.direction = Vector({1,1,0});
    camera.cs.basis.get_basis_coordinates(camera.direction);
    vector<vector<Entities::Ray>> Q = camera.get_rays_matrix(20, 80); // 50, 200
    Entities::Game::Canvas console(camera, 20,80);
    console.update(HE);
    //console.distances.print(); // trouble's here
    console.drawObjects(HE, 0.05f);
}