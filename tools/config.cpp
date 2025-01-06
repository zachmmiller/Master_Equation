#include "config.h"

lua_State* Lua_Init() {
    lua_State* L = luaL_newstate();
    luaL_openlibs(L);
    return L;
}

size_t Lua_Get_Vector_Size(lua_State* L) {
    if (lua_istable(L, -1)) {
        return lua_rawlen(L, -1);
    } else {
        return 0;
    }
}

void Lua_Pull_Vector_Entry(lua_State* L, int i) {
    lua_pushinteger(L, i);
    lua_gettable(L, -2);
}

void Lua_Pop_Vector_Entry(lua_State* L) { lua_pop(L, 1); }

template <typename T>
int Lua_Load_1d_Vector(lua_State* L, std::vector<T>& out) {
    out.clear();
    size_t d1 = Lua_Get_Vector_Size(L);
    if (d1 == 0) {
        return 1;
    } else {
        for (int i = 1; i < d1 + 1; i++) {
            Lua_Pull_Vector_Entry(L, i);
            if (lua_isnumber(L, -1)) {
                out.push_back((T)lua_tonumber(L, -1));
                Lua_Pop_Vector_Entry(L);
            } else {
                Lua_Pop_Vector_Entry(L);
                return 1;
            }
        }
        return 0;
    }
}

template <typename T>
int Lua_Load_Number(lua_State* L, T* out) {
    if (lua_isnumber(L, -1)) {
        *out = (T)lua_tonumber(L, -1);
        return 0;
    } else {
        return 1;
    }
}

int Lua_Load_Bool(lua_State* L, bool* out) {
    if (lua_isboolean(L, -1)) {
        *out = lua_toboolean(L, -1);
        return 0;
    } else {
        return 1;
    }
}

//
// Explicit template instantiations
// You must add to this list if you want additional templates
// For example, say you want to load a vector of shorts (16 bit integers), you would write:
//
// template int Lua_Load_1d_Vector<int_16t>(lua_State* L, std::vector <int_16t>& out);
//
// That's it! The compiler will figure out the rest.
// Note: int_16t and short are interchangeable on most architectures.
//

template int Lua_Load_1d_Vector<int>(lua_State* L, std::vector<int>& out);
template int Lua_Load_1d_Vector<double>(lua_State* L, std::vector<double>& out);
template int Lua_Load_Number<int>(lua_State* L, int* out);
template int Lua_Load_Number<double>(lua_State* L, double* out);
