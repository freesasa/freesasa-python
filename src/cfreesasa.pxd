from libc.stdio cimport FILE

cdef extern from "freesasa.h":
    ctypedef enum freesasa_algorithm:
        FREESASA_LEE_RICHARDS, FREESASA_SHRAKE_RUPLEY

    ctypedef enum freesasa_verbosity:
        FREESASA_V_NORMAL, FREESASA_V_NOWARNINGS, FREESASA_V_SILENT, FREESASA_V_DEBUG

    ctypedef enum freesasa_atom_class:
        FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_UNKNOWN

    cdef const int FREESASA_SUCCESS
    cdef const int FREESASA_FAIL
    cdef const int FREESASA_WARN

    cdef const int FREESASA_INCLUDE_HETATM
    cdef const int FREESASA_INCLUDE_HYDROGEN
    cdef const int FREESASA_SEPARATE_CHAINS
    cdef const int FREESASA_SEPARATE_MODELS
    cdef const int FREESASA_JOIN_MODELS
    cdef const int FREESASA_HALT_AT_UNKNOWN
    cdef const int FREESASA_SKIP_UNKNOWN

    cdef const int FREESASA_MAX_SELECTION_NAME

    ctypedef struct freesasa_parameters:
        freesasa_algorithm alg
        double probe_radius
        int shrake_rupley_n_points
        int lee_richards_n_slices
        int n_threads

    ctypedef struct freesasa_result:
        double total
        double *sasa
        int n_atoms

    ctypedef struct freesasa_nodearea:
        const char *name
        double total
        double main_chain
        double side_chain
        double polar
        double apolar
        double unknown

    ctypedef struct freesasa_node:
        pass

    ctypedef enum freesasa_nodetype:
        pass

    ctypedef struct freesasa_classifier:
        pass

    ctypedef struct freesasa_structure:
        pass

    ctypedef struct freesasa_selection:
        pass

    cdef extern const freesasa_parameters freesasa_default_parameters
    cdef extern const freesasa_classifier freesasa_default_classifier
    cdef extern const freesasa_classifier freesasa_residue_classifier
    cdef extern const freesasa_classifier freesasa_naccess_classifier
    cdef extern const freesasa_classifier freesasa_protor_classifier
    cdef extern const freesasa_classifier freesasa_oons_classifier

    freesasa_result* freesasa_calc_structure(const freesasa_structure *structure,
                                             const freesasa_parameters *parameters)

    freesasa_result* freesasa_calc_coord(const double *xyz,
                                         const double *radii,
                                         int n,
                                         const freesasa_parameters *parameters)

    void freesasa_result_free(freesasa_result *result)

    freesasa_classifier* freesasa_classifier_from_file(FILE *file)

    void freesasa_classifier_free(freesasa_classifier *classifier)

    int freesasa_structure_chain_residues(const freesasa_structure *structure,
                                          char chain,
                                          int *first,
                                          int *last)

    double freesasa_classifier_radius(const freesasa_classifier *classifier,
                                      const char *res_name,
                                      const char *atom_name)

    freesasa_atom_class freesasa_classifier_class(const freesasa_classifier *classifier,
                                                  const char *res_name,
                                                  const char *atom_name)

    const char* freesasa_classifier_class2str(freesasa_atom_class the_class)

    freesasa_selection * freesasa_selection_new(const char *command,
                                                const freesasa_structure *structure,
                                                const freesasa_result *result)

    void freesasa_selection_free(freesasa_selection *selection)

    const char * freesasa_selection_name(const freesasa_selection* selection)

    const char * freesasa_selection_command(const freesasa_selection* selection)

    double freesasa_selection_area(const freesasa_selection* selection)

    int freesasa_selection_n_atoms(const freesasa_selection* selection)

    int freesasa_write_pdb(FILE *output,
                           freesasa_result *result,
                           const freesasa_structure *structure)

    int freesasa_per_residue_type(FILE *output,
                                  freesasa_result *result,
                                  const freesasa_structure *structure)

    int freesasa_per_residue(FILE *output,
                             freesasa_result *result,
                             const freesasa_structure *structure)

    int freesasa_set_verbosity(freesasa_verbosity v)

    freesasa_verbosity freesasa_get_verbosity()

    freesasa_structure* freesasa_structure_from_pdb(FILE *pdb,
                                                    const freesasa_classifier* classifier,
                                                    int options)

    freesasa_structure** freesasa_structure_array(FILE *pdb,
                                                  int *n,
                                                  const freesasa_classifier* classifier,
                                                  int options)

    freesasa_structure* freesasa_structure_get_chains(
                                                const freesasa_structure *structure,
                                                const char *chains)

    freesasa_structure* freesasa_structure_new()

    int freesasa_structure_n(freesasa_structure *structure)

    void freesasa_structure_free(freesasa_structure* structure)

    const double* freesasa_structure_radius(const freesasa_structure *structure)

    void freesasa_structure_set_radius(freesasa_structure *structure,
                                       const double *radii)


    int freesasa_structure_add_atom(freesasa_structure *structure,
                                    const char* atom_name,
                                    const char* residue_name,
                                    const char* residue_number,
                                    char chain_label,
                                    double x, double y, double z)

    int freesasa_structure_add_atom_wopt(freesasa_structure *structure,
                                         const char* atom_name,
                                         const char* residue_name,
                                         const char* residue_number,
                                         char chain_label,
                                         double x, double y, double z,
                                         const freesasa_classifier *classifier,
                                         int options)

    const char* freesasa_structure_atom_name(const freesasa_structure *structure,
                                             int i)

    const char* freesasa_structure_atom_res_name(const freesasa_structure *structure,
                                                 int i)

    const char* freesasa_structure_atom_res_number(const freesasa_structure *structure,
                                                   int i)

    double freesasa_structure_atom_radius(const freesasa_structure *structure,
                                          int i)

    void freesasa_structure_atom_set_radius(const freesasa_structure *structure,
                                            int i,
                                            double radius)

    char freesasa_structure_atom_chain(const freesasa_structure *structure, int i)

    const double* freesasa_structure_coord_array(const freesasa_structure *structure)


    freesasa_nodearea freesasa_result_classes(const freesasa_structure *structure,
                                              const freesasa_result *result)

    freesasa_node * freesasa_tree_init(const freesasa_result *result,
                                       const freesasa_structure *structure,
                                       const char *name)

    int freesasa_node_free(freesasa_node *root)

    const freesasa_nodearea * freesasa_node_area(const freesasa_node *node)

    freesasa_node * freesasa_node_children(freesasa_node *node)

    freesasa_node * freesasa_node_next(freesasa_node *node)

    freesasa_node * freesasa_node_parent(freesasa_node *node)

    freesasa_nodetype freesasa_node_type(const freesasa_node *node)

    const char * freesasa_node_name(const freesasa_node *node)

    const char * freesasa_node_classified_by(const freesasa_node *node)

    int freesasa_node_atom_is_polar(const freesasa_node *node)

    int freesasa_node_atom_is_mainchain(const freesasa_node *node)

    double freesasa_node_atom_radius(const freesasa_node *node)

    const char * freesasa_node_atom_pdb_line(const freesasa_node *node)

    const char * freesasa_node_residue_number(const freesasa_node *node)

    int freesasa_node_residue_n_atoms(const freesasa_node *node)

    const freesasa_nodearea * freesasa_node_residue_reference(const freesasa_node *node)

    int freesasa_node_chain_n_residues(const freesasa_node *node)

    int freesasa_node_structure_n_chains(const freesasa_node *node)

    int freesasa_node_structure_n_atoms(const freesasa_node *node)

    const char * freesasa_node_structure_chain_labels(const freesasa_node *node)

    int freesasa_node_structure_model(const freesasa_node *node)

    const freesasa_result * freesasa_node_structure_result(const freesasa_node *node)

cdef extern from "freesasa_internal.h":
    int freesasa_write_pdb(FILE *output, freesasa_node *structure)