#pragma once

#include <limits>
#include <functional>
#include <thread>
#include <algorithm>
#include <vector>

namespace element_validity {
    void par_for(const size_t size, const size_t n_threads, const std::function<void(int, int, int)> &func)
    {
        if (n_threads == 1)
            func(0, size, /*thread_id=*/0); // actually the full for loop
        else
        {
            std::vector<std::thread> threads(n_threads);

            for (size_t t = 0; t < n_threads; t++)
            {
                threads[t] = std::thread(std::bind(
                    func,
                    t * size / n_threads,
                    (t + 1) == n_threads ? size : (t + 1) * size / n_threads,
                    t));
            }
            std::for_each(threads.begin(), threads.end(), [](std::thread &x) { x.join(); });
        }
    }
}
