#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../GenomeRearrangementSRC/simulator.h"
#include "../GenomeRearrangementSRC/genomes.h"


namespace py = pybind11;

PYBIND11_MODULE(GenomeRearrangement, m) {
    m.doc() = R"pbdoc(
        Genome Rearrangement simulator
        -----------------------

        .. currentmodule:: GenomeRearrangement

        .. autosummary::
           :toctree: _generate

           init_sim
           simulate_based_on_tree
    )pbdoc";

    py::class_<Simulator>(m, "Sim")
        .def(py::init<const std::string&>())
        .def("init_sim", &Simulator::initSim)
        .def("run_sim", &Simulator::simulateBasedOnTree)
        .def("get_events_count", &Simulator::get_event_counter_vectors)
        .def("get_tree", &Simulator::getSimTree)
        .def("set_seed", &Simulator::setSeed);

    py::class_<genomes>(m, "genomes")
        .def(py::init<std::string, std::string>())
        .def(py::init<const vector<genomeType> &, tree>())
        .def("get_sum_stats", &genomes::get_summary_stats_vector);

    py::class_<tree>(m, "tree")
        .def(py::init<>());
}
