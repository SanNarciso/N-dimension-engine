#ifndef ENGINE_ENTITIES_H
#define ENGINE_ENTITIES_H
#define ASPECT_RATIO 16/9
#endif //ENGINE_ENTITIES_H
#include <set>
#include <string>
#include <map>
#include <any>
#include "../Math/base.h"
#include <random>
#include <unordered_map>
#include <functional>
#include <cmath>

namespace Entities{
    class Ray{
    public:
        CoordinateSystem coord_sys;
        Point initial_point;
        Vector direction;
        Ray(CoordinateSystem cs, Point initial_pt, Vector dir) {
            coord_sys = cs;
            cs.basis.get_basis_coordinates(dir);
            initial_point = initial_pt;
            direction = dir;
        }
    };


    template <typename T>
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

    template <typename T>
    std::set<T> Identifier<T>::identifiers;


    class Entity {
    public:
        CoordinateSystem cs;
        Identifier<int> identifier;
        std::unordered_map<std::string, any> properties;
        Entity(const CoordinateSystem& crs = CoordinateSystem()) {
            cs = crs;
        }

        void set_property(const std::string& prop, const any& value) {
            properties[prop] = value;
        }

        any get_property(const std::string& prop) const {
            if (properties.count(prop) > 0) {
                return properties.at(prop);
            } else {
                return any();
            }
        }

        void remove_property(const std::string& prop) {
            properties.erase(prop);
        }

        any& operator[](const std::string& prop) {
            return properties[prop];
        }
    };


    class EntitiesList {
    public:
        std::vector<Entity> entities;
        void append(const Entity& entity) {
            entities.push_back(entity);
        }

        void remove(const Identifier<int>& id) {
            entities.erase(
                    std::remove_if(entities.begin(), entities.end(),
                                   [&id](const Entity& entity) { return entity.identifier.value == id.value; }),
                    entities.end());
        }

        Entity* get(const Identifier<int>& id) {
            auto it = std::find_if(entities.begin(), entities.end(),
                                   [&id](const Entity& entity) { return entity.identifier.value == id.value; });

            if (it != entities.end()) {
                return &(*it);
            }

            return nullptr;
        }

        void exec(std::function<void(Entity&)> func) {
            for (Entity& entity : entities) {
                func(entity);
            }
        }

        Entity& operator[](const Identifier<int>& id) {
            Entity* entity = get(id);
            if (entity) {
                return *entity;
            }
            throw std::out_of_range("Entity not found.");
        }
    };


    class Game {
    public:
        CoordinateSystem cs;
        EntitiesList entities;
        Game(CoordinateSystem crs, EntitiesList eentities){
            cs = crs;
            entities = eentities;
        }

        void run(){}
        void update(){}
        void exit(){}


    class Object : public Entities::Entity {
        public:
            Point position;
            Vector direction;
            Object(Point pos = Point(), Vector dir = Vector()) {
                cs.basis.get_basis_coordinates(dir);
                Vector normalized_dir((dir/dir.length).coordinates);
                cs.basis.get_basis_coordinates(normalized_dir);
                (*this).set_property("Position", pos);
                (*this).set_property("Direction", normalized_dir);
                position = pos;
                direction = normalized_dir;
            }

            void move(Vector direction) {
                this->position = position + direction;
            }

            void planar_rotate(int axis_1, int axis_2, float angle) {
                this->direction = Matrix::planar_rotation(axis_1, axis_2, angle, angle)*direction;
                this->position = Point(((Matrix::planar_rotation(axis_1, axis_2, angle, angle)*Vector(position.coordinates)).transpose()).matrix[0]);
            }

            void rotate_3d(vector<int> angles) {
                this->direction = Matrix::rotation_x(angles[0])*Matrix::rotation_y(angles[1])*Matrix::rotation_z(angles[2])*direction;
                this->position = Point((((Matrix::rotation_x(angles[0])*Matrix::rotation_y(angles[1])*Matrix::rotation_z(angles[2]))*Vector(position.coordinates)).transpose()).matrix[0]);
            }

            void set_position(Point position) {
                this->position = position;
            }

            void set_direction(Vector direction) {
                this->direction = direction;
            }
        };

        class Camera : public Object {
        public:
            float fov;
            float vfov;
            float draw_distance;
            Point look_at;
            Camera(float fov = 60, float draw_distance = INT_MAX)
                    : fov(toRadians(fov)), vfov(calculateVFOV(fov)), draw_distance(draw_distance) {}

            Camera(float fov, float vfov, float draw_distance)
                    : fov(toRadians(fov)), vfov(toRadians(vfov)), draw_distance(draw_distance) {}

            Camera(float fov, Point look_at, float draw_distance)
                    : fov(toRadians(fov)), vfov(calculateVFOV(fov)), draw_distance(draw_distance), look_at(look_at) {}

            Camera(float fov, float vfov, Point look_at, float draw_distance)
                    : fov(toRadians(fov)), vfov(toRadians(vfov)), draw_distance(draw_distance), look_at(look_at) {}



        private:
            float toRadians(float degrees) {
                return degrees * M_PI / 180.0;
            }

            float calculateVFOV(float fov) {
                return atan(tan(fov/2)*ASPECT_RATIO);
            }
        };

        template <typename EntityType>
        EntityType* get_entity_class() {
            return new EntityType(cs);
        }

        template <typename RayType>
        RayType* get_ray_class() {
            return new RayType(cs);
        }
    };


}