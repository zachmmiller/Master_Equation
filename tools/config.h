#ifndef CONFIG_H
#define CONFIG_H

#include <filesystem>
#include <string>
#include <vector>

#include "../vendor/Lua/src/lua.hpp"

size_t Lua_Get_Vector_Size(lua_State* L);
void Lua_Pull_Vector_Entry(lua_State* L, int i);
void Lua_Pop_Vector_Entry(lua_State* L);

template <typename T>
int Lua_Load_1d_Vector(lua_State* L, std::vector<T>& out);

template <typename T>
int Lua_Load_Number(lua_State* L, T* out);

int Lua_Load_Bool(lua_State* L, bool* out);

int Lua_Load_String(lua_State* L, std::string* out);

#endif
