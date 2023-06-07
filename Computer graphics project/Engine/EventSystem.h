#ifndef ENGINE_EVENTSYSTEM_H
#define ENGINE_EVENTSYSTEM_H

#endif //ENGINE_EVENTSYSTEM_H
#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <functional>
#include <any>

class EventSystem{
    std::unordered_map<std::string, vector<std::function<std::any(vector<std::any>)>>> events;

    void add_and_handle(const std::string name, const std::vector<function<any(vector<any>)>> handler){
        events[name] = handler;
    }

    void remove_and_handle(const std::string name){
        events.erase(name);
    }

    void trigger(std::string name, std::vector<any> arguments){
        for (auto func: events[name]){
            func(arguments);
        }
    }

    auto get_hadnler(std::string name){
        return events[name];
    }

    auto operator[](string name){
        return events[name];
    }
};