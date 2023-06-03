#ifndef ENGINE_ENTITIES_H
#define ENGINE_ENTITIES_H
#define ASPECT_RATIO 16/9
#endif //ENGINE_ENTITIES_H
#pragma once
#include <set>
#include <string>
#include "../../ncurses"
#include <map>
#include <any>
#include "../Math/base.h"
#include "../../Engine/EventSystem.h"
#include <random>
#include <unordered_map>
#include <functional>
#include <cmath>

namespace Entities {


    class Ray {
    public:
        CoordinateSystem coord_sys;
        Point initial_point;
        Vector direction;

        Ray(CoordinateSystem cs = CoordinateSystem(), Point initial_pt = Point(), Vector dir = Vector()) {
            coord_sys = cs;
            cs.basis.get_basis_coordinates(dir);
            initial_point = initial_pt;
            direction = dir;
        }

        void normalize() {
            Vector normalized = direction / direction.length;
            direction = normalized;
        }
    };


    template<typename T>
    class Identifier {
    public:
        static std::set<T> identifiers;
        T value;

        Identifier() {
            value = __generate__();
            identifiers.insert(value);
        }

        [[nodiscard]] T get_value() const {
            return value;
        }

    private:

        static T __generate__() {
            // Генерация нового значения идентификатора
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<T> dist(std::numeric_limits<T>::min(), std::numeric_limits<T>::max());

            T new_value = dist(gen);
            while (identifiers.count(new_value) > 0) {
                new_value = dist(gen);
            }

            return new_value;
        }
    };

    template<typename T>
    std::set<T> Identifier<T>::identifiers;


    class Entity {
    public:
        CoordinateSystem cs;
        Identifier<int> identifier;
        std::unordered_map<std::string, any> properties;

        Entity(const CoordinateSystem &crs = CoordinateSystem()) {
            cs = crs;
        }

        virtual double intersection_distance(Ray r){return 0;}

        void set_property(const std::string &prop, const any &value) {
            properties[prop] = value;
        }

        any get_property(const std::string &prop) const {
            if (properties.count(prop) > 0) {
                return properties.at(prop);
            } else {
                return any();
            }
        }

        void remove_property(const std::string &prop) {
            properties.erase(prop);
        }

        any &operator[](const std::string &prop) {
            return properties[prop];
        }
    };


    class EntitiesList {
    public:
        std::vector<Entity> entities;

        void append(const Entity &entity) {
            entities.push_back(entity);
        }

        void remove(const Identifier<int> &id) {
            entities.erase(
                    std::remove_if(entities.begin(), entities.end(),
                                   [&id](const Entity &entity) { return entity.identifier.value == id.value; }),
                    entities.end());
        }

        Entity *get(const Identifier<int> &id) {
            auto it = std::find_if(entities.begin(), entities.end(),
                                   [&id](const Entity &entity) { return entity.identifier.value == id.value; });

            if (it != entities.end()) {
                return &(*it);
            }

            return nullptr;
        }

        void exec(std::function<void(Entity &)> func) {
            for (Entity &entity : entities) {
                func(entity);
            }
        }

        Entity &operator[](const Identifier<int> &id) {
            Entity *entity = get(id);
            if (entity) {
                return *entity;
            }
            throw std::out_of_range("Entity not found.");
        }
    };


    EntitiesList Entities_List;


    class Game {
    public:
        CoordinateSystem cs;
        EntitiesList entities;
        EventSystem eventSystem;
        Game(CoordinateSystem crs, EntitiesList eentities) {
            cs = crs;
            Game::entities = eentities;
        }

        Game(CoordinateSystem crs, EntitiesList eentities, EventSystem es){
            cs = crs;
            Game::entities = eentities;
            eventSystem = es;
        }

        void run() {}

        void update() {}

        void exit() {}


        EventSystem get_event_system(){
            return (*this).eventSystem;
        }

        void apply_configuration();

        class Object : public Entities::Entity {
        public:
            Point position;
            Vector direction;

            Object(Point pos = Point(), Vector dir = Vector()) {
                cs.basis.get_basis_coordinates(dir);
                Vector normalized_dir((dir / dir.length).coordinates);
                cs.basis.get_basis_coordinates(normalized_dir);
                (*this).set_property("Position", pos);
                (*this).set_property("Direction", normalized_dir);
                position = pos;
                direction = normalized_dir;
            }

            void move(Vector dir) {
                this->position = position + dir;
            }

            void planar_rotate(int axis_1, int axis_2, float angle) {
                this->direction = Matrix::planar_rotation(axis_1, axis_2, angle, angle) * direction;
                cs.basis.get_basis_coordinates(this->direction);
                this->position = Point(((Matrix::planar_rotation(axis_1, axis_2, angle, angle) *
                                         Vector(position.coordinates)).transpose()).matrix[0]);
            }

            void rotate_3d(vector<int> angles) {
                this->direction =
                        Matrix::rotation_x(angles[0]) * Matrix::rotation_y(angles[1]) * Matrix::rotation_z(angles[2]) *direction;
                cs.basis.get_basis_coordinates(this->direction);
                this->position = Point((((Matrix::rotation_x(angles[0]) * Matrix::rotation_y(angles[1]) *
                                          Matrix::rotation_z(angles[2])) *
                                         Vector(position.coordinates)).transpose()).matrix[0]);
            }

            void set_position(Point pos) {
                this->position = pos;
            }

            void set_direction(Vector dir) {
                this->direction = dir;
                cs.basis.get_basis_coordinates(this->direction);
            }
        };



        class HyperPlane : public Game::Object {
        public:
            Point position;
            Vector normal;
            HyperPlane(Point pos, Vector norm) {
                cs.basis.get_basis_coordinates(norm);
                (*this).set_property("Position", pos);
                (*this).set_property("Normal vector", norm);
                position = pos;
                normal = norm;
                cs.basis.get_basis_coordinates(normal);
                Entities_List.append((*this));
            }

            void planar_rotate(int axis_1, int axis_2, float angle){
                Vector res(((Matrix::planar_rotation(axis_1, axis_2, angle, angle)*normal).transpose()).matrix[0]);
                cs.basis.get_basis_coordinates(res);
                normal = res;
                cs.basis.get_basis_coordinates(normal);
                (*this).set_property("Normal vector", normal);
            }

            void rotate_3d(double alpha, double beta, double gamma){
                Vector res((((Matrix::rotation_x(alpha) * Matrix::rotation_y(beta) * Matrix::rotation_z(gamma)) * normal).transpose()).matrix[0]);
                cs.basis.get_basis_coordinates(res);
                normal = res;
                cs.basis.get_basis_coordinates(normal);
                (*this).set_property("Normal vector", normal);
            }

            double intersection_distance(Ray ray){
                cs.basis.get_basis_coordinates(ray.direction);
                // if direction is parallel to the plane, it never intersects
                if ((normal.transpose()*cs.basis.Gramm*ray.direction).matrix[0][0] == 0){
                    return -1;
                }
                double t_hat = -((normal.transpose()*cs.basis.Gramm*(Vector(ray.initial_point.coordinates) - Vector((*this).position.coordinates))).matrix[0][0])/(normal.transpose()*cs.basis.Gramm*ray.direction).matrix[0][0];
                Vector intersection = Vector(ray.initial_point.coordinates) + ray.direction*t_hat;
                Vector diff(intersection.dim);
                for (int i = 0; i < intersection.dim; ++i){
                    diff.coordinates[i] = intersection.coordinates[i] - ray.initial_point.coordinates[i];
                }
                cs.basis.get_basis_coordinates(diff);
                return diff.length;
            }
        };


        class HyperEllipsoid : public Game::Object {
        public:
            Point position; Vector direction; vector<Vector> semiaxes;

            HyperEllipsoid(Point pos, Vector dir, vector<Vector> ss){
                (*this).set_property("Position", pos);
                (*this).set_property("Direction", dir);
                (*this).set_property("Semiaxes", semiaxes);
                cs.basis.get_basis_coordinates(dir);
                position = pos; direction = dir; semiaxes = ss;
                for (auto &q : semiaxes){
                    cs.basis.get_basis_coordinates(q);
                }
                Entities_List.append((*this));
            }

            void rotate_3d(double alpha, double beta, double gamma){
                Vector new_direction((((Matrix::rotation_x(alpha)*Matrix::rotation_y(beta)*Matrix::rotation_z(gamma))*direction).transpose()).matrix[0]);
                cs.basis.get_basis_coordinates(new_direction);
                (*this).set_property("Direction", new_direction);
                for (int i = 0; i < semiaxes.size(); ++i){
                    semiaxes[i] = Vector(((((Matrix::rotation_x(alpha)*Matrix::rotation_y(beta)*Matrix::rotation_z(gamma))*semiaxes[i]).transpose()).matrix[0]));
                }
                (*this).set_property("Semiaxes", semiaxes);
            }

            double intersection_distance(Ray r){
                double A = 0;
                double B = 0;
                double C = 0;
                cs.basis.get_basis_coordinates(r.direction);
                for (int i = 0; i < cs.basis.dimension_of_vector_space; ++i){
                    A += pow(r.initial_point.coordinates[i]/semiaxes[i].coordinates[i],2);              // here we just expect that ellipsoid is in its canon basis
                    B += (2*r.direction.coordinates[i]*r.initial_point.coordinates[i]) / pow(semiaxes[i].coordinates[i],2);
                    C += pow(r.initial_point.coordinates[i] / semiaxes[i].coordinates[i], 2);
                }
                if ((B - 4*A*(C-1)) < 0){
                    return -1;
                }
                double found_t = min((-B + sqrt(B - 4*A*(C-1)))/(2*A), (-B - sqrt(B - 4*A*(C-1)))/(2*A));
                Vector intersection_point((r.initial_point + r.direction*found_t).coordinates);
                Vector diff(intersection_point.dim);
                for (int i = 0; i < diff.dim; ++i){
                    diff.coordinates[i] = intersection_point.coordinates[i] - r.initial_point.coordinates[i];
                }
                cs.basis.get_basis_coordinates(diff);
                return diff.length;
            }
        };



        class Camera : public Object {
        public:
            // Vector direction is not needed here as
            float fov;
            float vfov;
            float draw_distance;
            optional<Point> look_at;

            Camera(float fov = 60, float draw_distance = INT_MAX)
                    : fov(toRadians(fov)), vfov(calculateVFOV(fov)), draw_distance(draw_distance) {Entities_List.append((*this));}

            Camera(float fov, float vfov, float draw_distance)
                    : fov(toRadians(fov)), vfov(toRadians(vfov)), draw_distance(draw_distance) {Entities_List.append((*this));}

            Camera(float fov, Point look_at, float draw_distance)
                    : fov(toRadians(fov)), vfov(calculateVFOV(fov)), draw_distance(draw_distance), look_at(look_at) {Entities_List.append((*this));}

            Camera(float fov, float vfov, Point look_at, float draw_distance)
                    : fov(toRadians(fov)), vfov(toRadians(vfov)), draw_distance(draw_distance), look_at(look_at) {Entities_List.append((*this));}


            vector<vector<Ray>> get_rays_matrix(int n, int m) {

                vector<vector<Ray>> A(n, vector<Ray>(m));
                if (look_at.has_value()) {
                    for (int i = 0; i < look_at->coordinates.size(); ++i) {
                        direction.coordinates[i] = look_at->coordinates[i] - position.coordinates[i];
                    }
                }
                double alpha = fov;
                double beta = vfov;
                double d_alpha = alpha / n;
                double d_beta = beta / m;
                vector<double> alpha_i(n, 0);
                vector<double> beta_j(m, 0);

                for (int i = 0; i < n; ++i) {
                    alpha_i[i] = d_alpha * i - alpha / 2;
                }
                for (int j = 0; j < m; ++j) {
                    beta_j[j] = d_beta * j - beta / 2;
                }

                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < m; ++j) {
                        Vector v_i_j((((Matrix::rotation_z(alpha_i[i]) * Matrix::rotation_y(beta_j[j])) *
                                       direction).transpose()).matrix[0]);
                        cs.basis.get_basis_coordinates(v_i_j);
                        Vector _v_i_j = v_i_j * (pow(direction.length, 2) /
                                                 ((direction.transpose() * cs.basis.Gramm * v_i_j).matrix[0][0]));
                        //A[i][j].direction = _v_i_j;
                        cs.basis.get_basis_coordinates(_v_i_j);
                        A[i][j] = Ray(cs, position, _v_i_j);
                    }
                }
                return A;
            }

        private:
            float toRadians(float degrees) {
                return degrees * M_PI / 180.0;
            }

            float calculateVFOV(float fov) {
                return atan(tan(fov / 2) * ASPECT_RATIO);
            }
        };

        class Canvas{
        public:
            int n;
            int m;
            Matrix distances = Matrix(1920, 1080);

            Canvas(int nn = 1920, int mm = 1080) {
                n = nn;
                m = mm;
                distances = Matrix(n,m);
            }

            void draw(){}/// NCURSES, DRAWING EVERYTHING

            void update(Camera cam){
                /*
                 * what is needed to be done here?
                 * first, we send rays from camera and thus get matrix with rays_from_camera
                 * then in each cell of the gained matrix we need to find intersection of the ray contained
                 * in the current matrix cell and any Entity from Entites list, of no intersection found, set the cell eqaul to draw_distance
                 */
                vector<vector<Ray>> camera_ray_cast = cam.get_rays_matrix(n, m);
                for (int i = 0; i < n; ++i){
                    for (int j = 0; j < m; ++j){
                        for (int e = 0; e < Entities_List.entities.size(); ++e){
                            double dist_ray_i_j = Entities_List.entities[e].intersection_distance(camera_ray_cast[i][j]);
                            if (dist_ray_i_j != -1){
                                distances.matrix[i][j] = dist_ray_i_j;
                            } else {
                                distances.matrix[i][j] = cam.draw_distance;
                            }
                        }
                    }
                }
            }
        };


        class Console: public Canvas{
        public:
            string charmap = ".:;><+r*zsvfwqkP694VOGbUAKXH8RD#$B0MNWQ%&@";
            void drawObjects() {
                initscr(); // Initialize ncurses
                raw(); // Disable line buffering
                keypad(stdscr, TRUE); // Enable special key handling
                noecho(); // Don't echo input

                int rows, cols;
                getmaxyx(stdscr, rows, cols); // Get terminal size

                double maxDistance = 0.0;
                for (const auto& row : canvas.distances) {
                    for (double distance : row) {
                        if (distance > maxDistance) {
                            maxDistance = distance;
                        }
                    }
                }

                for (int i = 0; i < rows; ++i) {
                    for (int j = 0; j < cols; ++j) {
                        double distance = canvas.distances[i][j];
                        char symbol;

                        if (distance <= 0.0) {
                            symbol = '#'; // Object is close
                        } else {
                            double normalizedDistance = distance / maxDistance;
                            if (normalizedDistance < 0.25) {
                                symbol = ' ';
                            } else if (normalizedDistance < 0.5) {
                                symbol = '.';
                            } else if (normalizedDistance < 0.75) {
                                symbol = '-';
                            } else {
                                symbol = '=';
                            }
                        }

                        mvaddch(i, j, symbol);
                    }
                }

                refresh(); // Refresh the screen
                getch(); // Wait for user input
                endwin(); // End ncurses
            }
        };

//        class Configuration{
//
//        };

        template<typename EntityType>
        EntityType *get_entity_class() {
            return new EntityType(cs);
        }

        template<typename RayType>
        RayType *get_ray_class() {
            return new RayType(cs);
        }
    };
}

