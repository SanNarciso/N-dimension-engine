#ifndef ENGINE_EVENTSYSTEM_H
#define ENGINE_EVENTSYSTEM_H

#endif //ENGINE_EVENTSYSTEM_H
#pragma once
#include <vector>
#include <string>

class EventSystem{
    std::unordered_map<std::string, vector<std::function<any(vector<any>)>>> events;

    void add_and_handle(const string name, const vector<function<any(vector<any>)>> handler){
        events[name] = handler;
    }

    void remove_and_handle(const string name){
        events.erase(name);
    }

    void trigger(string name, vector<any> arguments){
        for (auto func: events[name]){
            func(arguments);
        }
    }

    auto get_hadnler(string name){
        return events[name];
    }

    auto operator[](string name){
        return events[name];
    }
};