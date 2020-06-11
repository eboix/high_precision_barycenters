#include "skgeom.hpp"
#include "funcs.hpp"

#include <CGAL/Rotational_sweep_visibility_2.h>
#include <CGAL/Arr_segment_traits_2.h>

#include <CGAL/Arrangement_2.h>

////#include <CGAL/CORE_algebraic_number_traits.h>
////#include <CGAL/Arr_Bezier_curve_traits_2.h>

typedef Segment_Arrangement_2::Ccb_halfedge_circulator          Ccb_halfedge_circulator;
typedef Segment_Arrangement_2::Isolated_vertex_iterator         Isolated_vertex_iterator;

typedef PySequenceCirculator<Ccb_halfedge_circulator> PyCcbHalfedgeCirculator;


Halfedge_handle insert_non_intersecting_curve_in_arr(Segment_Arrangement_2& arr, Segment_2& seg) {
    return CGAL::insert_non_intersecting_curve(arr, seg);
}

void insert_non_intersecting_curves_in_arr(Segment_Arrangement_2& arr, std::vector<Segment_2>& segs) {
    CGAL::insert_non_intersecting_curves(arr, segs.begin(), segs.end());
}

// HACK: method added by Enric.
void insert_segments_in_arr(Segment_Arrangement_2& arr, std::vector<Segment_2>& segs) {
    CGAL::insert(arr, segs.begin(), segs.end());
}

Vertex_handle insert_point_in_arr(Segment_Arrangement_2 & arr, Point_2& p) {
    Vertex_handle h = CGAL::insert_point(arr, p);
    return h;
}

bool remove_vertex_from_arr(Segment_Arrangement_2& arr, Vertex_handle v) {
    return CGAL::remove_vertex(arr, v);
}

Face_handle remove_edge_from_arr(Segment_Arrangement_2& arr, Halfedge_handle h) {
    return CGAL::remove_edge(arr, h);
}

py::object find_in_arrangement(Segment_Arrangement_2& arr, Point_2& query_point) {
    CGAL::Arr_naive_point_location<Segment_Arrangement_2> pl(arr);
    CGAL::Arr_point_location_result<Segment_Arrangement_2>::Type obj = pl.locate(query_point);

    const Vertex_const_handle*   v;
    const Halfedge_const_handle* e;
    const Face_const_handle*     f;

    // Casting to non - const types since they are only defined in skgeom
    if ((f = boost::get<Face_const_handle>(&obj)))
    {
        // Note it's very important to return `Handles` here.
        return py::cast(*(Face_handle*) f);
    }
    else if ((e = boost::get<Halfedge_const_handle>(&obj)))
    {
        return py::cast(*(Halfedge_handle*) e);   
    }
    else if ((v = boost::get<Vertex_const_handle>(&obj)))
    {
        return py::cast(*(Vertex_handle*) v);
    }
    return py::none();
}

namespace skgeom
{
    using namespace pybind11;

    template <typename Iterator, return_value_policy Policy>
    struct circulator_state {
        Iterator it;
        bool first;
    };

    template <return_value_policy Policy = return_value_policy::reference_internal,
              typename HandleType = void,
              typename Iterator,
              typename Sentinel,
              typename ValueType = decltype(*std::declval<Iterator>()),
              typename... Extra>
    iterator make_handle_iterator(Iterator first, Sentinel last, Extra &&... extra) {
        typedef detail::iterator_state<Iterator, Sentinel, false, Policy> state;

        if (!detail::get_type_info(typeid(state), false)) {
            class_<state>(handle(), "iterator", pybind11::module_local())
                .def("__iter__", [](state &s) -> state& { return s; })
                .def("__next__", [](state &s) -> HandleType {
                    if (!s.first_or_done)
                        ++s.it;
                    else
                        s.first_or_done = false;
                    if (s.it == s.end) {
                        s.first_or_done = true;
                        throw stop_iteration();
                    }
                    // do not dereference here, just cast to handle
                    return HandleType(s.it);
                }, std::forward<Extra>(extra)..., Policy);
        }

        return cast(state{first, last, true});
    }

    template <return_value_policy Policy = return_value_policy::reference_internal,
              typename HandleType = void,
              typename Iterator,
              typename ValueType = decltype(*std::declval<Iterator>()),
              typename... Extra>
    iterator make_circulator(Iterator first, Extra &&... extra) {
        typedef circulator_state<Iterator, Policy> state;

        if (!detail::get_type_info(typeid(state), false)) {
            class_<state>(handle(), "circulator", pybind11::module_local())
                .def("__iter__", [](state &s) -> state& { return s; })
                .def("__next__", [](state &s) -> HandleType {
                    if (s.first) { return HandleType(s.it); }
                    else { return HandleType(++s.it); };
                }, std::forward<Extra>(extra)..., Policy);
        }

        return cast(state{first});
    }

    template <return_value_policy Policy = return_value_policy::reference_internal,
              typename Iterator,
              typename Sentinel,
              typename... Extra>
    iterator make_hole_iterator(Iterator first, Sentinel last, Extra &&... extra) {
        typedef detail::iterator_state<Iterator, Sentinel, false, Policy> state;

        if (!detail::get_type_info(typeid(state), false)) {
            class_<state>(handle(), "iterator", pybind11::module_local())
                .def("__iter__", [](state &s) -> state& { return s; })
                .def("__next__", [](state &s) {
                    if (!s.first_or_done)
                        ++s.it;
                    else
                        s.first_or_done = false;
                    if (s.it == s.end) {
                        s.first_or_done = true;
                        throw stop_iteration();
                    }
                    return skgeom::make_circulator<Policy, Halfedge_handle>(*s.it);
                }, std::forward<Extra>(extra)..., Policy);
        }

        return cast(state{first, last, true});
    }
}

void init_arrangement(py::module &m) {
    py::module sub = m.def_submodule("arrangement");

    py::class_<Segment_Arrangement_2>(sub, "Arrangement")
        .def(py::init<>())
        .def_property_readonly("halfedges", [](Segment_Arrangement_2& s) {
           return skgeom::make_handle_iterator<py::return_value_policy::reference_internal, Halfedge_handle>(s.halfedges_begin(), s.halfedges_end());
        }, py::keep_alive<0, 1>())
        .def_property_readonly("faces", [](Segment_Arrangement_2& s) {
           return skgeom::make_handle_iterator<py::return_value_policy::reference_internal, Face_handle>(s.faces_begin(), s.faces_end());
        }, py::keep_alive<0, 1>())
        .def_property_readonly("vertices", [](Segment_Arrangement_2& s) {
           return skgeom::make_handle_iterator<py::return_value_policy::reference_internal, Vertex_handle>(s.vertices_begin(), s.vertices_end());
        }, py::keep_alive<0, 1>())
        .def("insert_non_intersecting_curve", &insert_non_intersecting_curve_in_arr)
        .def("insert_non_intersecting_curves", &insert_non_intersecting_curves_in_arr)
        .def("insert_segments", &insert_segments_in_arr) // HACK: method added by Enric.
        .def("insert", [](Segment_Arrangement_2& self, Segment_2& segment) {
            CGAL::insert(self, segment);
        })
        .def("unbounded_face", static_cast<Face_handle (Segment_Arrangement_2::*)()>(&Segment_Arrangement_2::unbounded_face))
        // .def("insert_from_left_vertex", &Segment_Arrangement_2::insert_from_left_vertex)
        // .def("insert_from_right_vertex", &Segment_Arrangement_2::insert_from_right_vertex)
        // .def("insert_in_face_interior", &Segment_Arrangement_2::insert_in_face_interior)
        // .def("insert_at_vertices", &Segment_Arrangement_2::insert_at_vertices)
        // .def("split_edge", &Segment_Arrangement_2::split_edge)
        .def("address", [](Segment_Arrangement_2& self) { std::cout << (void*) &self << std::endl; })
        .def("insert_point", &insert_point_in_arr)
        .def("remove_vertex", &remove_vertex_from_arr)
        .def("remove_edge", &remove_edge_from_arr)
        .def("find", &find_in_arrangement)
    ;

    py::class_<Vertex, Vertex_handle>(sub, "Vertex")
        .def("point", [](Vertex& self) { return self.point(); })
        .def_property_readonly("incident_halfedges", [](Vertex& self) {
            return skgeom::make_circulator<py::return_value_policy::reference_internal, Halfedge_handle>(self.incident_halfedges());
        }, py::keep_alive<0, 1>())
    ;

    py::class_<Halfedge, Halfedge_handle>(sub, "Halfedge")
        .def("prev", static_cast<Halfedge_handle (Halfedge::*)()>(&Halfedge::prev))
        .def("next", static_cast<Halfedge_handle (Halfedge::*)()>(&Halfedge::next))
        .def("twin", static_cast<Halfedge_handle (Halfedge::*)()>(&Halfedge::twin))
        .def("source", static_cast<Vertex_handle (Halfedge::*)()>(&Halfedge::source))
        .def("target", static_cast<Vertex_handle (Halfedge::*)()>(&Halfedge::target))
        .def("curve", [](Halfedge& he) -> Segment_2 {
            return he.curve();
        })
        .def("face", static_cast<Face_handle (Halfedge::*)()>(&Halfedge::face))
        // .def_property_readonly("ccb", [](Halfedge& self) {
        //     return py::make_circulator(self.ccb());
        // }, py::keep_alive<0, 1>())
    ;

    py::class_<Face, Face_handle>(sub, "Face")
        .def("is_unbounded", [](Face & f) -> bool { return f.is_unbounded(); })
        .def("has_outer_ccb", &Face::has_outer_ccb)
        .def("number_of_holes", &Face::number_of_holes)
        .def("number_of_isolated_vertices", [](Face & f) -> int { return f.number_of_isolated_vertices(); })
        .def_property_readonly("isolated_vertices", [](Face& self) {
            return py::make_iterator(self.isolated_vertices_begin(), self.isolated_vertices_end()); 
        }, py::keep_alive<0, 1>())
        .def_property_readonly("outer_ccb", [](Face& self) {
            return skgeom::make_circulator<py::return_value_policy::reference_internal, Halfedge_handle>(self.outer_ccb()); 
        }, py::keep_alive<0, 1>())
        .def_property_readonly("holes", [](Face& self) {
            return skgeom::make_hole_iterator(self.holes_begin(), self.holes_end()); 
        }, py::keep_alive<0, 1>())
    ;

}
