from cfreesasa cimport *

class ResidueArea:
    """
    Stores absolute and relative areas for a residue

    Attributes:
        residueType (str): Type of Residue
        residueNumber (str): Residue number
        hasRelativeAreas (bool): False if there was noe reference area to calculate relative areas from

        total (float): Total SASA of residue
        polar (float): Polar SASA
        apolar (float): Apolar SASA
        mainChain (float): Main chain SASA
        sideChain (float): Side chain SASA

        relativeTotal (float): Relative total SASA
        relativePolar (float): Relative polar SASA
        relativeApolar (float): Relative Apolar SASA
        relativeMainChain (float): Relative main chain SASA
        relativeSideChain (float): Relative side chain SASA
    """

    residueType = ""
    residueNumber = ""
    hasRelativeAreas = False

    total = 0
    polar = 0
    apolar = 0
    mainChain = 0
    sideChain = 0

    relativeTotal = 0
    relativePolar = 0
    relativeApolar = 0
    relativeMainChain = 0
    relativeSideChain = 0


cdef class Result:
    """
    Stores results from SASA calculation.

    The type of object returned by :py:func:`freesasa.calc()`,
    not intended to be used outside of that context.
    """

    cdef freesasa_result* _c_result
    cdef freesasa_node* _c_root_node

    ## The constructor
    def __init__ (self):
        self._c_result = NULL
        self._c_root_node = NULL

    ## The destructor
    def __dealloc__(self):
        if self._c_result is not NULL:
            freesasa_result_free(self._c_result)
        if self._c_root_node is not NULL:
            freesasa_node_free(self._c_root_node)

    def nAtoms(self):
        """
        Number of atoms in the results.

        Returns:
            int: Number of atoms.
        """
        if self._c_result is not NULL:
            return self._c_result.n_atoms
        return 0

    def totalArea(self):
        """
        Total SASA.

        Returns:
            The total area in Å^2.
        Raises:
            AssertionError: If no results have been associated with the object.
        """
        assert(self._c_result is not NULL)
        return self._c_result.total

    def atomArea(self,i):
        """
        SASA for a given atom.

        Args:
            i (int): index of atom.

        Returns:
            float: SASA of atom i in Å^2.

        Raise:
            AssertionError: If no results have been associated
                      with the object or if index is out of bounds
        """
        assert(self._c_result is not NULL)
        assert(i < self._c_result.n_atoms)
        return self._c_result.sasa[i]

    def residueAreas(self):
        """
        Get SASA for all residues including relative areas if available for the
        classifier used.

        Returns dictionary of results where first dimension is chain label and
        the second dimension residue number. I.e. ``result["A"]["5"]`` gives the
        :py:class:`.ResidueArea` of residue number 5 in chain A.

        Relative areas are normalized to 1, but can be > 1 for
        residues in unusual conformations or at the ends of chains.

        Returns:
            dictionary

        Raise:
            AssertionError: If no results or structure has been associated
                 with the object.
        """
        assert(self._c_result is not NULL)
        assert(self._c_root_node is not NULL, "Result.residueAreas can only be called on results generated directly or indirectly by freesasa.calc()")

        cdef freesasa_node* result_node = <freesasa_node*> freesasa_node_children(self._c_root_node)
        cdef freesasa_node* structure = <freesasa_node*> freesasa_node_children(result_node)
        cdef freesasa_node* chain
        cdef freesasa_node* residue
        cdef freesasa_nodearea* c_area
        cdef freesasa_nodearea* c_ref_area

        result = {}

        chain = <freesasa_node*> freesasa_node_children(structure)
        while (chain != NULL):
            residue  = <freesasa_node*> freesasa_node_children(chain)
            chainLabel = freesasa_node_name(chain)
            result[chainLabel] = {}

            while (residue != NULL):
                c_area = <freesasa_nodearea*> freesasa_node_area(residue)
                c_ref_area = <freesasa_nodearea*> freesasa_node_residue_reference(residue)
                residueNumber = freesasa_node_residue_number(residue).strip()
                residueType = freesasa_node_name(residue).strip()

                area = ResidueArea()

                area.residueType = residueType
                area.residueNumber = residueNumber

                area.total = c_area.total
                area.mainChain = c_area.main_chain
                area.sideChain = c_area.side_chain
                area.polar = c_area.polar
                area.apolar = c_area.apolar

                if (c_ref_area is not NULL):
                    area.hasRelativeAreas = True
                    area.relativeTotal = self._safe_div(c_area.total, c_ref_area.total)
                    area.relativeMainChain = self._safe_div(c_area.main_chain, c_ref_area.main_chain)
                    area.relativeSideChain = self._safe_div(c_area.side_chain, c_ref_area.side_chain)
                    area.relativePolar = self._safe_div(c_area.polar, c_ref_area.polar)
                    area.relativeApolar = self._safe_div(c_area.apolar, c_ref_area.apolar)

                result[chainLabel][residueNumber] = area

                residue = <freesasa_node*> freesasa_node_next(residue)

            chain = <freesasa_node*> freesasa_node_next(chain)

        return result

    def write_pdb(self, filename):
        if self._c_root_node is NULL:
            raise AssertionError('Result root node points to NULL. Unable to write to a pdb file.')

        cdef freesasa_node *result_node    = <freesasa_node*> freesasa_node_children(self._c_root_node)
        cdef freesasa_node *structure_node = <freesasa_node*> freesasa_node_children(result_node)
        cdef freesasa_node *chain_node     = <freesasa_node*> freesasa_node_children(structure_node)
        cdef freesasa_node *residue_node   = <freesasa_node*> freesasa_node_children(chain_node)
        cdef freesasa_node *atom_node      = <freesasa_node*> freesasa_node_children(residue_node)

        cdef const char * atom_pdb_line    =  freesasa_node_atom_pdb_line(atom_node)

        if atom_pdb_line is NULL:
            raise AssertionError(
                "Atom PDB Line is NULL. You are probably trying to write a PDB file from a Bio.Structure."
            )
        
        cdef FILE *f = NULL

        f = fopen(filename, 'w')
        freesasa_write_pdb(f, self._c_root_node)
        fclose(f)

    def _safe_div(self,a,b):
        try:
            return a/b
        except ZeroDivisionError:
            return float('nan')

    def _get_address(self, size_t ptr2ptr):
        cdef freesasa_result **p = <freesasa_result**> ptr2ptr
        p[0] = self._c_result
