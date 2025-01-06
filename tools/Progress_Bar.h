//
// Created by Zach Miller on 12/1/23.
//

#ifndef PROGRESS_BAR_H
#define PROGRESS_BAR_H

#include <iostream>
#include <ostream>
#include <string>

class Progress_Bar {
   public:
    Progress_Bar(int N, int size, const char* name) : m_N(N), m_pos(-1), m_size(size), m_name(name), m_bracket("|"), m_todo(" "), m_done("█"){};
    ~Progress_Bar() = default;

    void Start() { Update(); }

    void Update() {
        if (m_pos >= m_N) {
            return;
        }

        m_pos++;
        int p = std::floor(100.0 * m_pos / m_N);
        int n_done_chars = std::floor(m_size * m_pos / m_N);
        int n_todo_chars = m_size - n_done_chars;
        std::cout << "\r" << m_bracket;

#if !defined(_WIN64) || !defined(_WIN32)
        for (int _ = 0; _ < n_done_chars; _++) {
            std::cout << m_done;
        }

        for (int _ = 0; _ < n_todo_chars; _++) {
            std::cout << m_todo;
        }
        std::cout << m_bracket;
#endif
        if (p < 10) {
            std::cout << "   " << std::to_string(p) << "%";
        } else if (p < 100) {
            std::cout << "  " << std::to_string(p) << "%";
        } else  // AKA p >= 100
        {
            std::cout << " " << std::to_string(p) << "%";
        }

        std::cout << " " << m_name << std::flush;
    }

    void End() { std::cout << std::endl; }

   private:
    int m_N;
    int m_pos;
    int m_size;

    std::string m_name;

    char m_bracket[2];
    char m_todo[2];
    char m_done[4];
};

#endif  // CDMS_ANALYSIS_PROGRESS_BAR_H
