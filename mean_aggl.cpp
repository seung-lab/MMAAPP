//
// Copyright (C) 2013-present  Aleksandar Zlateski <zlateski@mit.edu>
// ------------------------------------------------------------------
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <boost/heap/binomial_heap.hpp>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <vector>
#include <zi/disjoint_sets/disjoint_sets.hpp>

template <class T>
struct edge_t
{
    uint64_t v0, v1;
    T        w;
};

template <class T>
using region_graph = std::vector<edge_t<T>>;

template <class T, class C = std::greater<T>>
struct heapable_edge;

template <class T, class C = std::greater<T>>
struct heapable_edge_compare
{
    bool operator()(heapable_edge<T, C>* const a,
                    heapable_edge<T, C>* const b) const
    {
        C c;
        return c(b->edge.w, a->edge.w);
    }
};

template <class T, class C = std::greater<T>>
using heap_type = boost::heap::binomial_heap<
    heapable_edge<T, C>*, boost::heap::compare<heapable_edge_compare<T, C>>>;

template <class T, class C>
struct heapable_edge
{
    edge_t<T> edge;
    typename heap_type<T, C>::handle_type handle;
};

template <class T, class Compare = std::greater<T>, class Plus = std::plus<T>,
          class Limits = std::numeric_limits<T>>
inline std::vector<uint64_t> agglomerate(std::vector<edge_t<T>> const& rg,
                                         T const& threshold, uint64_t const n)
{
    Compare comp;
    Plus    plus;
    heap_type<T, Compare> heap;

    zi::disjoint_sets<uint64_t> sets(n);
    std::vector<std::map<uint64_t, heapable_edge<T, Compare>*>> incident(n);

    std::vector<heapable_edge<T, Compare>> edges(rg.size());
    for (std::size_t i = 0; i < rg.size(); ++i)
    {
        edges[i].edge                = rg[i];
        edges[i].handle              = heap.push(&edges[i]);
        incident[rg[i].v0][rg[i].v1] = &edges[i];
        incident[rg[i].v1][rg[i].v0] = &edges[i];
    }

    while (heap.size() && comp(heap.top()->edge.w, threshold))
    {
        auto e = heap.top();
        heap.pop();

        auto v0 = e->edge.v0;
        auto v1 = e->edge.v1;

        if (v0 != v1)
        {

            {
                auto s0 = sets.find_set(v0);
                auto s1 = sets.find_set(v1);
                auto s  = sets.join(s0, s1);

                // std::cout << "Joined " << s0 << " and " << s1 << " to " << s
                //           << " at " << e->edge.w << "\n";
                if (s0 != s1) {
                    std::cout << s0 << " " << s1 << " " << s << " " << e->edge.w << std::endl;
                }
            }

            if (incident[v0].size() > incident[v1].size())
            {
                std::swap(v0, v1);
            }

            // v0 is dissapearing from the graph

            // earase the edge e = {v0,v1}
            incident[v0].erase(v1);
            incident[v1].erase(v0);

            // loop over other edges e0 = {v0,v}
            for (auto& e0 : incident[v0])
            {
                incident[e0.first].erase(v0);
                if (incident[v1].count(e0.first)) // {v0,v} and {v1,v} exist, we
                                                  // need to merge them
                {
                    auto& e1   = incident[v1][e0.first]; // edge {v1,v}
                    e1->edge.w = plus(e1->edge.w, e0.second->edge.w);
                    heap.update(e1->handle);
                    {
                        // std::cout
                        //     << "Removing: " <<
                        //     sets.find_set(e0.second->edge.v0)
                        //     << " to " << sets.find_set(e0.second->edge.v1)
                        //     << " of " << e0.second->edge.w << " ";
                        e0.second->edge.w = Limits::max();
                        heap.increase(e0.second->handle);
                        heap.pop();
                    }
                }
                else
                {
                    if (e0.second->edge.v0 == v0)
                        e0.second->edge.v0 = v1;
                    if (e0.second->edge.v1 == v0)
                        e0.second->edge.v1 = v1;
                    incident[e0.first][v1] = e0.second;
                    incident[v1][e0.first] = e0.second;
                }
                incident[v0].erase(e0.first);
            }
        }
    }

    std::vector<uint64_t> remaps(n, std::numeric_limits<uint64_t>::max());

    uint64_t next = 0;

    for (uint64_t i = 0; i < n; ++i)
    {
        auto s = sets.find_set(i);
        if (remaps[s] == std::numeric_limits<uint64_t>::max())
        {
            remaps[s] = next++;
        }
        remaps[i] = remaps[s];
    }

    std::cout << "Total of " << next << " segments\n";

    while (heap.size())
    {
        auto e = heap.top();
        heap.pop();
        auto v0 = e->edge.v0;
        auto v1 = e->edge.v1;
        auto s0 = sets.find_set(v0);
        auto s1 = sets.find_set(v1);
        std::cout << s0 << " " << s1 << " " << e->edge.w << std::endl;
    }

    return remaps;
}

typedef struct atomic_edge
{
    uint64_t u1;
    uint64_t u2;
    uint64_t area;
    explicit constexpr atomic_edge(uint64_t w1 = 0, uint64_t w2 = 0, uint64_t a = 0)
        : u1(w1)
        , u2(w2)
        , area(a)
    {
    }
} atomic_edge_t;

struct mean_edge
{
    double sum;
    double num;
    atomic_edge_t * repr;

    explicit constexpr mean_edge(double s = 0, double n = 1, atomic_edge_t * r = NULL)
        : sum(s)
        , num(n)
        , repr(r)
    {
    }
};

struct mean_edge_plus
{
    mean_edge operator()(mean_edge const& a, mean_edge const& b) const
    {
        atomic_edge_t * new_repr = NULL;
        if (a.repr->area > b.repr->area)
            new_repr = a.repr;
        else
            new_repr = b.repr;
        return mean_edge(a.sum + b.sum, a.num + b.num, new_repr);
    }
};

struct mean_edge_greater
{
    bool operator()(mean_edge const& a, mean_edge const& b) const
    {
        return a.sum / a.num > b.sum / b.num;
    }
};

struct mean_edge_limits
{
    static constexpr mean_edge max()
    {
        return mean_edge(std::numeric_limits<double>::max(), 1, NULL);
    }
};

template <class CharT, class Traits>
::std::basic_ostream<CharT, Traits>&
operator<<(::std::basic_ostream<CharT, Traits>& os, mean_edge const& v)
{
    os << v.sum / v.num << " " << v.repr->u1 << " " <<  v.repr->u2;
    return os;
}

int main()
{
    std::vector<edge_t<mean_edge>> rg;

    std::size_t v, n;
    std::cin >> v >> n;

    for (std::size_t i = 0; i < n; ++i)
    {
        edge_t<mean_edge> e;
        std::cin >> e.v0 >> e.v1 >> e.w.sum >> e.w.num;
        atomic_edge_t * ae = new atomic_edge_t(e.v0, e.v1, e.w.sum);
        e.w.repr = ae;
        rg.push_back(e);
    }

    auto res = agglomerate<mean_edge, mean_edge_greater, mean_edge_plus,
                           mean_edge_limits>(rg, mean_edge(66, 1), v);

    for (auto& e : res)
    {
        std::cout << e << "\n";
    }
}
