import IMP
import IMP.core
import IMP.atom
import IMP.threading
import IMP.algebra
import IMP.pmi
import IMP.pmi.tools
import random
import numpy
import math
import IMP.pmi.restraints


class Loop():
    def __init__(self, start=0, length=0, se_before_id=-1, se_after_id=-1):
        self.start=start
        self.length=length
        self.se_before_id=se_before_id
        self.se_after_id=se_after_id
    def get_loop_residues(self):
        return list(range(self.start, self.start+self.length))

    def generate_loop_id(self):
        return str(self.start)+"_"+str(self.length)+"_"+str(self.se_before_id)+"_"+str(self.se_after_id)

    def is_res_in_loop(self, resnum):
        if resnum in self.get_loop_residues():
            return True
        else:
            return False

class SSEThread(IMP.ModelObject):
    def __init__(self, root_hier = None, sequences={}, max_res_dist=4.0, se_clash_dist=4.0):
        
        self.sequences = sequences # list of fastas
        
        if root_hier:
            self.root_hier = root_hier
            self.model = self.root_hier.get_model()
            self.structure_chain = IMP.atom.Chain.setup_particle(IMP.Particle(self.model), "X")
            #self.sequence_hierarchy = IMP.atom.Hierarchy.setup_particle(IMP.Particle(self.model))
            self._update_SSE_hier()
        else:
            self.model = IMP.Model()
            self._setup_SSE_system()
           
        
        self.max_res_dist = max_res_dist
        self.se_clash_dist = se_clash_dist
        
        if len(self.sequences) > 0:
            self.build_sequence_chains()

    def _setup_SSE_system(self):
        
        self.structure_chain = IMP.atom.Chain.setup_particle(IMP.Particle(self.model), "X")
        self.sequence_hierarchy = IMP.atom.Hierarchy.setup_particle(IMP.Particle(self.model))

        self.root_hier.add_child(self.structure_chain)
        self.root_hier.add_child(self.sequence_hierarchy)

    def _update_SSE_hier(self):
        self.sequence_hierarchy = []
        mols = IMP.atom.get_by_type(self.root_hier, IMP.atom.MOLECULE_TYPE)
        mol_names = {mol.get_name():i for i, mol in enumerate(mols)}
        
        for prot, seq in self.sequences.items():
            seq_chain = self.build_sequence_chains(seq = self.sequences[prot], prot = prot)
            if prot in mol_names.keys():
                self.root_hier.add_child(seq_chain)
                #mols[k].get_parent().add_child(seq_chain)
            else:
                self.root_hier.add_child(seq_chain)
            # Fix this
            self.sequence_hierarchy.append(seq_chain)
                    
        #for mol in mols:
        #    prot = mol.get_name()
        #    # Check of coordinates for protein to be threaded are already present
        #    if prot in self.sequences.keys():
        #        seq_chain = self.build_sequence_chains(seq = self.sequences[prot], prot = prot)
        #        #seq_chain = IMP.atom.Chain.setup_particle(IMP.Particle(self.model, "seq_chain_"+mol.get_name()), "S"+mol.get_name())
        #        mol.get_parent().add_child(seq_chain)
        #    else:
                
        
        #print(IMP.atom.show_with_representations(self.root_hier))
        mols = IMP.atom.get_by_type(self.root_hier, IMP.atom.MOLECULE_TYPE)
        
        
    def build_sequence_chains(self, seq = [], prot='1'):
        # Given a fasta sequence and the root hierarchy, build the sequence chain

        sse_res_id = 1

        if '.' in prot:
            chain_num = prot.split('.')[1]
            prot = prot.split('.')[0]
        else:
            chain_num = 0

        new_seq_chain = IMP.atom.Chain.setup_particle(IMP.Particle(self.model, "seq_chain_"+prot), "S"+prot)
        IMP.atom.Copy.setup_particle(new_seq_chain,int(chain_num))
        
        for r in range(len(seq)):
            pr = IMP.Particle(self.model)
            
            # Setup residue/mass/xyzr
            res = IMP.atom.Residue.setup_particle(pr,
                                                  IMP.pmi.tools.get_residue_type_from_one_letter_code(seq[r]),
                                                  sse_res_id)
            #res_particles.append(res.get_particle())
            IMP.atom.Mass.setup_particle(res.get_particle(), IMP.atom.get_mass(res.get_residue_type()))
            IMP.core.XYZR.setup_particle(res.get_particle())

            # Initialize coordinates to 0,0,0
            IMP.core.XYZ(pr).set_coordinates((0,0,0))
            
            # We use the "are optimized" flag to indicate structured or not
            # Here, we begin with an unstructured residue
            IMP.core.XYZ(pr).set_coordinates_are_optimized(False)
            
            new_seq_chain.add_child(res)
            sse_res_id+=1
            # Add to hierarchy
            #self.add_child(new_seq_chain)

        #self.compute_chain_breaks()
        return new_seq_chain
                
    
    def get_sequence_hierarchy(self):
        return self.sequence_hierarchy

    def get_updated_hierarchy(self):
        return self.root_hier

        
    def compute_chain_breaks(self):
        self.chain_breaks=[]
        start_res=0
        for chain in self.sequence_hierarchy.get_children():
            self.chain_breaks.append(len(chain.get_children())+start_res)
            start_res = self.chain_breaks[-1]
        return self.chain_breaks

    def generate_structure_elements_from_hier(self, hier, pdb = '', elements = [], sec_struct = [], stride=None):
        '''
        Given a hierarchy and a list of SE, create the structure elements
        '''


        if len(elements) != len(sec_struct):
            print(elements, sec_struct)
            raise ValueError("Number of elements and secondary structure assignments are different")
            
        self.structure_elements = {}
        for i, e in enumerate(elements):
            m_temp = IMP.Model()
            h = IMP.atom.read_pdb(pdb, m_temp, IMP.atom.CAlphaPDBSelector())
            cas = IMP.atom.Selection(h,
                                     chain_id=e[0],
                                     residue_indexes=range(e[1],e[1]+e[2]),
                                     atom_type=IMP.atom.AT_CA).get_selected_particles()
            coords = []
            for p in cas:
                coords.append(IMP.core.XYZ(p).get_coordinates())
                
            if stride is not None:
                self.structure_elements[i] = self.setup_structure_element(coords, sec_struct[i], 0, 1, 0, stride)
            else:
                self.structure_elements[i] = self.setup_structure_element(coords, sec_struct[i], 0, 1, 0, stride = None)

        del m_temp

        return self.structure_elements

    def add_coordinates_to_sequence_chain(self, hier, pdb, prot, copy_index = 0, fixed=[]):
        for i, f in enumerate(fixed):
            m_temp = IMP.Model()
            h = IMP.atom.read_pdb(pdb, m_temp, IMP.atom.CAlphaPDBSelector())
            for r in range(f[1],f[1]+f[2]):    
                p_pdb = IMP.atom.Selection(h,
                                           chain_id=f[0],
                                           residue_index=r,
                                           atom_type=IMP.atom.AT_CA).get_selected_particles()
            
                p_seq = IMP.atom.Selection(hier,
                                           chain_id='S'+prot,
                                           residue_index=r,
                                           copy_index = i).get_selected_particles()
                if len(p_pdb)==1 and len(p_seq)==1:
                    IMP.core.XYZ(p_seq[0]).set_coordinates(IMP.core.XYZ(p_pdb[0]).get_coordinates())
                    IMP.core.XYZ(p_seq[0]).set_coordinates_are_optimized(True)
                    
        del m_temp
        

    def extract_structure_elements_from_smotif_pdbs(self, smotif_pdbs, stride=None):
        '''
        Given a set of PDBs, create two structure elements for each
        '''
        m2 = IMP.Model()
        self.structure_elements = {}
        se_ix = 0
        for p in smotif_pdbs:
            h = IMP.atom.read_pdb(p, m2, IMP.atom.CAlphaPDBSelector())

            # First, get CA atoms
            cas = IMP.atom.Selection(h, atom_type=IMP.atom.AT_CA).get_selected_particles()

            coords0 = []
            coords1 = []
            first = True

            residm1 = 0

            # First number of residues is in filename
            # A_1_26_43_44_56_0.pdb
            fname = p.split("/")[-1]
            nres0 = int(fname.split("_")[3])-int(fname.split("_")[2]) + 1
            #nres1 = int(p.split("_")[5])-int(p.split("_")[4]) + 1

            # Collect coordinates
            for c in range(nres0):
                coords0.append(IMP.core.XYZ(cas[c]).get_coordinates())
            for d in range(c+1, len(cas)):
                coords1.append(IMP.core.XYZ(cas[d]).get_coordinates())

            if stride is not None:
                n1 = int(fname.split("_")[2])
                n2 = int(fname.split("_")[3])+1
                n3 = int(fname.split("_")[4])
                n4 = int(fname.split("_")[5])+1
                stride0 = [stride[r] for r in range(n1, n2)]
                stride1 = [stride[r] for r in range(n3, n4)]
                #print(stride0, stride1)
                self.structure_elements[se_ix] = self.setup_structure_element(coords0, "H", 0, 1, 0, stride0)
                self.structure_elements[se_ix+1] = self.setup_structure_element(coords1, "H", 0, 1, 0, stride1)
            else:
                self.structure_elements[se_ix] = self.setup_structure_element(coords0, "H", 0, 1, 0)
                self.structure_elements[se_ix+1] = self.setup_structure_element(coords1, "H", 0, 1, 0)

            #print(p)
            #print("   ", structure_elements[se_ix].get_coordinates()[0:3])
            #print("   ", structure_elements[se_ix+1].get_coordinates()[0:3])
            se_ix+=2

        #print("Keys", self.structure_elements.keys())
        return self.structure_elements

    def setup_structure_element(self, coordinates, sec_struct, start_residue=0, polarity=1, offset=0, chain='0', stride=None):

        se_pi = IMP.Particle(self.model)
        se_hier = IMP.atom.Hierarchy.setup_particle(se_pi)
        self.root_hier.add_child(se_hier)

        # Add coordinates to hierarchy
        for c in range(len(coordinates)):
            coord = coordinates[c]
            np = IMP.Particle(self.model)
            hp = IMP.atom.Hierarchy.setup_particle(np)
            xyz = IMP.core.XYZR.setup_particle(np)
            xyz.set_coordinates(coord)
            xyz.set_radius(1.0)                     # radius doesn't really matter atm.
            IMP.atom.Mass.setup_particle(np, 1.0)   # Mass doesn't matter either
            se_hier.add_child(hp)

            if stride is not None:
                IMP.atom.SecondaryStructureResidue.setup_particle(np, stride[c][0], stride[c][1], stride[c][2])

        se = IMP.threading.StructureElement.setup_particle(self.model, se_pi.get_index(),
                                                           start_residue,
                                                           polarity,
                                                           len(coordinates),
                                                           offset,
                                                           chain)

        if sec_struct=="H":
            IMP.atom.SecondaryStructureResidue.setup_particle(se_pi, 1, 0, 0)
        elif sec_struct=="S":
            IMP.atom.SecondaryStructureResidue.setup_particle(se_pi, 0, 1, 0)
        elif sec_struct=="C":
            IMP.atom.SecondaryStructureResidue.setup_particle(se_pi, 0, 0, 1)
        else:
            raise Exception("Secondary structure designation must be H, S or C.", sec_struct, "was provided")
        return se

    def construct_se_min_loop_table(self):
        # Now build the hash table
        sids = list(self.structure_elements.keys())
        self.se_min_loop_table = numpy.zeros((len(sids), len(sids)))

        for sid0 in range(len(sids)):
            si0 = self.structure_elements[sids[sid0]]
            for sid1 in range(sid0, len(sids)):
                si1 = self.structure_elements[sids[sid1]]

                if self.are_structure_elements_clashing(si0, si1):
                    val = -1
                else:
                    val = self.get_minimum_spanning_residues(si0.get_coordinates()[-1], si1.get_coordinates()[0])
                    self.se_min_loop_table[sid0][sid1] = int(val)
                    
                    val = self.get_minimum_spanning_residues(si1.get_coordinates()[-1], si0.get_coordinates()[0])
                    self.se_min_loop_table[sid1][sid0] = int(val)

        return self.se_min_loop_table 

    def get_minimum_spanning_residues(self, xyz1, xyz2):
        return math.ceil(IMP.algebra.get_distance(xyz1, xyz2) / self.max_res_dist)

    def are_structure_elements_clashing(self, se0, se1):
        dists = []
        for c0 in se0.get_coordinates():
            for c1 in se1.get_coordinates():
                dist = IMP.algebra.get_distance(c0, c1)
                dists.append(dist)
                if dist < self.se_clash_dist:
                    return True
        return False

    def get_built_residues(self):
        built_res = []
        for p in IMP.atom.Selection(self.sequence_hierarchy).get_selected_particles():
            if IMP.core.XYZ(p).get_coordinates_are_optimized():
                built_res.append(p)
        return built_res

    def get_built_structure_element_ids(self, sort=False):
        built_ids = []
        built_sr = []

        # Use numpy.where instead
        for s in self.structure_elements.keys():
            sr = self.start_res_list[s]#int(self.structure_elements[s].get_start_res())
            if sr != 0: 
                built_ids.append(s)

        if sort:
            return self.sort_seids(built_ids)
        else:
            return built_ids

    def get_start_res_list(self):
        self.start_res_list=[]
        for se in range(len(self.structure_elements.keys())):
            self.start_res_list.append(int(self.structure_elements[se].get_start_res()))
        return self.start_res_list

    def get_all_loops(self):
        return self.get_loops(self.get_built_structure_element_ids())

    def sort_seids(self, seids):
        if seids is None or len(seids)==0:
            return []

        ses = []
        
        for s in seids:
            ses.append((s, self.structure_elements[s].get_start_res()))

        sorted_ids = sorted(ses, key=lambda x: x[1])
        return [s[0] for s in sorted_ids]

    def get_loops(self, seids, set_as=True):
        # Given the current state of the system, return a list of loops

        sorted_ses_ids = self.sort_seids(seids)

        loops = {}

        new_loop_start = 1
        chain = 0

        new_se_before_id = -1
        for s in range(len(sorted_ses_ids)):
            seid = sorted_ses_ids[s]
            start_res = int(self.structure_elements[seid].get_start_res())

            # if this start_res is after a chain, then we need to 
            if int(start_res) >= self.chain_breaks[chain]:
                start = new_loop_start
                length = self.chain_breaks[chain]-new_loop_start-1
                new_loop = Loop(start=start,
                                length=length,
                                se_before_id=new_se_before_id,
                                se_after_id=-1)
                loop_id = new_loop.generate_loop_id()
                loops[loop_id] = new_loop
                            
                chain+=1

                # New loop startes at N-term of the new chain
                new_loop_start=self.chain_breaks[chain]
                new_se_before_id = -1
                s-=1 # Need to do this SE again

            else:
                new_loop = Loop(start=new_loop_start,
                                length=start_res-new_loop_start,
                                se_before_id=new_se_before_id,
                                se_after_id=seid)
                loop_id = new_loop.generate_loop_id()
                loops[loop_id] = new_loop
                new_loop_start = start_res+self.loops[loop_id].length
                new_se_before_id = seid

        # Now, manually do the last loop(s)
        for i in range(chain, len(self.chain_breaks)):
            length = self.chain_breaks[i]-new_loop_start+1
            new_loop = Loop(start=new_loop_start,
                        length=length,
                        se_before_id=new_se_before_id,
                        se_after_id=-1)
            loop_id = new_loop.generate_loop_id()
            loops[loop_id] = new_loop
            new_loop_start=self.chain_breaks[i]+1
            new_se_before_id=-1

        if set_as:
            self.loops=loops

        return loops

    def is_se_allowed_in_model(self, structure_element_id):
        built_seids = self.get_built_structure_element_ids()

        # If structure element is built, then it is allowed
        if structure_element_id in built_seids:
            return True

        # If any built SE clashes (table=-1), then it is not allowed
        for b in built_seids:
            if self.se_min_loop_table[structure_element_id][b]==-1:
                return False
        return True

    def get_available_start_residues(self, structure_element_id, exclude_self=False, same_loop=False):
        # Given a structure element, find the start residues available to it
        available_residues = []

        se = self.structure_elements[structure_element_id]
        # If not allowed due to clash, then return nothing
        if not self.is_se_allowed_in_model(structure_element_id):
            return []

        # If allowed, then look at all loops
        se_len = int(se.get_length())

        if exclude_self:
            built_ids = self.get_built_structure_element_ids()
            if structure_element_id in built_ids:
                loops = self.get_loops(self.get_built_structure_element_ids().remove(structure_element_id), set_as=False)
            else:
                loops = self.loops
        else:
            loops = self.loops

        for lk in loops.keys():

            loop = loops[lk]

            # If we only want residues in the same loop, return them here
            if same_loop:
                if se.get_start_res() > loop.start and se.get_start_res() < loop.start+loop.length:
                    br = int(self.se_min_loop_table[loop.se_before_id][structure_element_id])
                    ar = int(self.se_min_loop_table[structure_element_id][loop.se_before_id])
                    resis = loop.get_loop_residues()
                    return resis[br:len(resis)-ar-loop.length]
                else:
                    continue

            # If loop is not large enough to hold SE, then continue
            if loop.length < se_len:
                continue

            # if loop starts at N-terminus, then there is no loop before offset
            if loop.se_before_id==-1:
                before_res=0
            # Otherwise, see how many residues this SE needs to be from the SE before
            else:
                before_res = int(self.se_min_loop_table[loop.se_before_id][structure_element_id])
            # Now do the same for the C-terminus
            if loop.se_after_id==-1:
                after_res=0
            else:
                after_res = int(self.se_min_loop_table[structure_element_id][loop.se_before_id])

            # If, after subtracting the offsets, the loop is now too small, continue without adding residues
            if loop.length -  before_res - after_res < se_len:
                continue

            else:
                available_residues+=list(range(loop.start+before_res, loop.start+loop.length-after_res-se_len))

        return available_residues

    def get_loop_ids_bracketing_seid(self, seid):
        # From the loops dictionary, return the loop before and loop after the given SEID
        # If there is no loop before or after (i.e., it's at a terminus), return None for that element

        before_loop_id = None
        after_loop_id = None

        for lk in self.loops.keys():
            loop = self.loops[lk]
            if loop.se_before_id == seid:
                after_loop_id = lk
            if loop.se_after_id == seid:
                before_loop_id = lk
        return before_loop_id, after_loop_id

    def get_residue_chain_termini(self, resnum):
        # For a given residue, return the N- and C-terminal residues for the chain it is in
        chains = [0]+self.chain_breaks
        for c in range(1, len(chains)):
            if resnum < chains[c]:
                return [chains[c-1]+1, chains[c]]

        raise Exception("Residue", resnum, "not in sequence")

    def remove_se_from_loop_table(self, seid):
        
        # First, get the two loops that we will combine        
        before_loop_id, after_loop_id = self.get_loop_ids_bracketing_seid(seid)
      
        #print("---", [i for i in self.loops.keys()])
        #print("---", self.structure_elements[seid].get_start_res(), before_loop_id, after_loop_id)
        
        # If no loops are found, this SEID was not built into the model
        # Just return.
        if before_loop_id is None and after_loop_id is None:
            return

        # Compute parameters for the new loop
        if before_loop_id is not None:
            before_loop = self.loops[before_loop_id]
            start = before_loop.start
            se_before_id = before_loop.se_before_id
            #print("BLI:", before_loop_id)
        else:
            # This SE starts at the N-terminus of a chain
            # Start of the new loop will now be the start residue of the SEID
            # And there is no SE before this loop now (since it starts at N-terminus
            after_loop = self.loops[after_loop_id]
            start = self.get_residue_chain_termini(after_loop.start)[0]
            se_before_id = -1

            #print(" - start:", start, after_loop.start, self.get_residue_chain_termini(after_loop.start))
        if after_loop_id is not None:
            after_loop = self.loops[after_loop_id]
            length = after_loop.start + after_loop.length - start
            #print("**", length, after_loop.start, after_loop.length, start)
            se_after_id = after_loop.se_after_id
        else:
            # This SE ends at the C-terminus.
            #print("C-term")
            length = self.get_residue_chain_termini(start)[1] - start
            se_after_id = -1
        new_loop = Loop(start=start,
                        length=length,
                        se_before_id=se_before_id,
                        se_after_id=se_after_id)

        new_loop_id = new_loop.generate_loop_id()

        #print("   REMOVE",before_loop_id, after_loop_id, new_loop_id)

        # Remove the two loops we connected
        if before_loop_id is not None:
            del self.loops[before_loop_id]
        if after_loop_id is not None:
            del self.loops[after_loop_id]

        # And add our new loop to the dictionary
        self.loops[new_loop_id] = new_loop

    def add_se_to_loop_table(self, seid, new_sr):
        # If the new SR is zero, do not add anything
        if new_sr == 0:
            return

        this_loop = None
        # Find the loop containing new_sr
        #print("--LOOPS--")
        for lk in self.loops.keys():
            #print(" ",lk, new_sr, "|",self.loops[lk].start,self.loops[lk].length,"|", len(self.loops[lk].get_loop_residues()), self.loops[lk].is_res_in_loop(new_sr))
            #print(" ",new_sr, self.loops[lk].start, self.loops[lk].get_loop_residues()[-1], self.loops[lk].is_res_in_loop(new_sr))
            if self.loops[lk].is_res_in_loop(new_sr):
                this_loop = self.loops[lk]
                break
        
        if this_loop is None:
            raise Exception("No loop found for SE starting at", new_sr, [(self.loops[s].start, self.loops[s].start+self.loops[s].length) for s in self.loops.keys()])
        mode = None
        # if the new start res is the start of this loop, then we only create one new loop
        if new_sr == this_loop.start:
            start = new_sr+self.structure_elements[seid].get_length()
            new_loop = Loop(start=start,
                            length=this_loop.start+this_loop.length-start,
                            se_before_id=seid,
                            se_after_id=this_loop.se_after_id)
            self.loops[new_loop.generate_loop_id()] = new_loop
            mode = 0
        # If the new SE goes all the way to the end of this loop, then again, only create one loop
        elif new_sr == this_loop.start+this_loop.length-self.structure_elements[seid].get_length():
            start = this_loop.start
            length = new_sr-start
            se_before_id = this_loop.se_before_id
            se_after_id = seid
            new_loop = Loop(start=start,
                           length=length,
                           se_before_id=se_before_id,
                           se_after_id=se_after_id)
            self.loops[new_loop.generate_loop_id()] = new_loop
            se_before_id = this_loop.se_before_id
            mode = 1
        # Otherwise, we need to create two loops out of one.
        
        else:
            # Make the first loop
            start = this_loop.start
            length = new_sr - start
            se_before_id = this_loop.se_before_id
            se_after_id = seid
            new_loop = Loop(start=start,
                           length=length,
                           se_before_id=se_before_id,
                           se_after_id=se_after_id)
            self.loops[new_loop.generate_loop_id()] = new_loop
            start = new_sr + self.structure_elements[seid].get_length()
            length = this_loop.start+this_loop.length-start
            se_before_id = seid
            se_after_id = this_loop.se_after_id

            new_loop = Loop(start=start,
                           length=length,
                           se_before_id=se_before_id,
                           se_after_id=se_after_id)
            self.loops[new_loop.generate_loop_id()] = new_loop
            mode = 2
        if mode == 2:
            snd_loop = self.loops[list(self.loops.keys())[-2]]
            #print("-22- New Loops", new_loop.start, new_loop.length, "|", snd_loop.start, snd_loop.length, "|", new_sr, self.structure_elements[seid].get_length())
        
        else:
            i=0
            #print("---- New Loop", new_loop.start, new_loop.length, "|", new_sr, self.structure_elements[seid].get_length())
        del self.loops[lk]



    def modify_loop_table(self, seid, new_sr, remove=True):
        
        # First, remove the seid loop
        self.remove_se_from_loop_table(seid)

        # Second, add the new loops to the table
        self.add_se_to_loop_table(seid, new_sr)

    def create_se_connectivity_restraints(self, structure_elements=None, function=None, n_sds=2):
        # Make N-1 connectivity restraints.

        if structure_elements is None:
            structure_elements = self.structure_elements
        # All restraints will not be used at all times.  
        # Any SECR restraints not used will be assigned dummy particles
        # which indicates the restraint should be evaluated as zero

        # make dummy particles for when we want to zero out SECR restraints
        self.secr_dummy_p0 = IMP.Particle(self.model)
        self.secr_dummy_p1 = IMP.Particle(self.model)
        
        if function is None:
            function = IMP.core.HarmonicUpperBound(0, 0.2) 

        self.secrs = []
        for i in range(len(structure_elements)-1):
            secr = IMP.threading.StructureElementConnectivityRestraint(self.model,
                                                                       function,
                                                                       self.secr_dummy_p0.get_index(), 
                                                                       self.secr_dummy_p1.get_index(),
                                                                       n_sds)

            self.secrs.append(secr)
        return self.secrs

    def update_SECR_restraints(self):
        
        try:
            n_secr = len(self.secrs)
        #if there are no secr restraints, then just exit
        except:
            return

        secr_id = 0
    
        # Find the connected SEs. Easiest through loop table
        for lk in self.loops.keys():
            if self.loops[lk].se_before_id == -1 or self.loops[lk].se_after_id == -1:
                continue
            else:
                # Those loops with SEs before and after indicate a connection
                se0 = self.structure_elements[self.loops[lk].se_before_id]
                se1 = self.structure_elements[self.loops[lk].se_after_id]

                # Assign these seids to the restraint
                self.secrs[secr_id].assign_particles(se0.get_particle_index(), se1.get_particle_index())

                secr_id+=1
        
        # Place the dummy particles in the remaining SECRs
        for i in range(secr_id, n_secr):
            self.secrs[secr_id].assign_particles(self.secr_dummy_p0.get_index(), self.secr_dummy_p1.get_index())

    def update_system(self):
        self.update_SECR_restraints()

    def add_SS_to_sequence(self, psipred):
        self.sequence_ss = {}
        for res in psipred.keys():
            sys_part = IMP.atom.Selection(self.sequence_hierarchy, residue_index=res).get_selected_particle_indexes()[0]
            IMP.atom.SecondaryStructureResidue.setup_particle(self.model, sys_part, psipred[res][0], psipred[res][1], psipred[res][2])
         
            self.sequence_ss[res] = psipred[res]

    def add_SS_to_ses(self, secstruct):
        for seid in secstruct.keys():
            se_res = self.structure_elements[seid].get_children()
            for res in range(len(se_res)):
                ss = secstruct[seid][res]
                IMP.atom.SecondaryStructureResidue.setup_particle(self.model, se_res[res].get_particle_index(), 
                        ss[0], ss[1], ss[2])

    def get_all_model_ss_propensities(self):
        # Returns a dictionary of tuples
        self.model_ss = {}

        built_ses = self.get_built_structure_element_ids()
        all_resids = set(self.sequence_ss.keys())
        for seid in built_ses:
            se = self.structure_elements[seid]
            se_res = IMP.atom.Hierarchy(se).get_children()
            resids = range(se.get_start_res(), se.get_start_res()+se.get_length())
            # Remove these built resids
            [all_resids.remove(r) for r in resids]
            for r in resids:
                self.model_ss[r] = IMP.atom.SecondaryStructureResidue(self.model, 
                                    se_res[r-se.get_start_res()].get_particle_index()).get_all_probabilities()
        # All unbuilt residues are modeled as coil
        for r in all_resids:
            self.model_ss[r]=(0,0,1)

        return self.model_ss

class CompletenessRestraint(IMP.Restraint):
    '''
    A restraint on the number of residues in the model that are built

    Evaluated as the number of unbuilt residues time a constant (slope)

    Unbuilt residues are determined by X = 0.0, as is used in SSEThread
    '''
    def __init__(self, root_hierarchy, slope=1.0):
        self.model = root_hierarchy.get_model()
        self.xyzs = [IMP.core.XYZ(p) for p in IMP.atom.Selection(root_hierarchy, resolution=1).get_selected_particles()]
        self.slope = slope
        IMP.Restraint.__init__(self, self.model, "MockRestraint %1%")

    def get_number_of_built_residues(self):
        built=0
        for xyz in self.xyzs:
            if xyz.get_coordinate(0) != 0.0:
                built+=1

        return built

    def unprotected_evaluate(self, da):
        return ((len(self.xyzs)-self.get_number_of_built_residues())*self.slope)

    def do_get_inputs(self):
        return []

    def do_show(self, fh):
        fh.write('CompletenessRestraint')


class SecondaryStructureParsimonyRestraint(IMP.Restraint):
    def __init__(self, system, psipred):
        self.system = system
        #self.system.add_SS_to_ses(se_sec_struct)
        self.system.add_SS_to_sequence(psipred)
        IMP.Restraint.__init__(self, self.system.model, "SecondaryStructureParsimonyRestraint")

    def unprotected_evaluate(self, da):
        score = 0
        # First compute structured residues, as determined by built structure elements:
        model_ss = self.system.get_all_model_ss_propensities()
        
        for res in self.system.sequence_ss:
            resscore=0
            for i in range(3):
                resscore+=self.system.sequence_ss[res][i]*model_ss[res][i]
            if resscore < 0.1:
                resscore = 0.1
            score+=-1*math.log(resscore)

        return score

    def do_get_inputs(self):
        return []
    
    def do_show(self, fh):
        fh.write('SecondaryStructureParsimonyRestraint')



class SSEAdditionMover(IMP.threading.SSEThreadMover):
    def __init__(self, system, structure_elements):
        '''
        A move that proposes deletion of a random SE from the model
        @param system :: the SSEThread object
        @param structure_elements :: A dictionary of structure elements with key - SEID
        '''
        se_pis = [structure_elements[seid].get_particle_index() for seid in structure_elements.keys()]
        #se_pis = [se.get_particle_index() for se in structure_elements]
        IMP.threading.SSEThreadMover.__init__(self, 
                system.model, 
                se_pis,
                system.sequence_hierarchy.get_particle_index())
        self.structure_elements = structure_elements
        self.system = system

    def get_se(self, seid):
        return self.structure_elements[self.seid]

    def get_random_structure_element(self):
        self.structure_elements[numpy.randint(0,len(structure_elements))]

    def do_propose(self, new_start_res=None, seid=None):
        # Get a random structure element that is not built
        if seid is None:
            rand_srs = list(self.structure_elements.keys())
            random.shuffle(rand_srs)

        #print(rand_srs)
            for seid in rand_srs:
                se = self.structure_elements[seid]
                sr = se.get_start_res()
                #print("  ", seid, sr)
                if sr == 0:
                    break

                # If none are zero, we can't add one.  Retnurn nothing.
                if seid == rand_srs[-1]:
                    return IMP.core.MonteCarloMoverResult([se.get_particle_index()], 0.0)

             #print("Add", seid, sr)
        
        self.seid = seid # Store the seid in a persisten variable
        se = self.structure_elements[self.seid]

        # Pick a random available start residue and change the SE start_res to that
        if new_start_res is None:
            available_start_res = self.system.get_available_start_residues(seid)
            if len(available_start_res)==0:
                return IMP.core.MonteCarloMoverResult([se.get_particle_index()], 0.0)

            new_start_res = random.choice(available_start_res)
        se.set_start_res_key(new_start_res)
        
        # Update the model
        se_pix = self.structure_elements[self.seid].get_particle_index()
        self.transform_coordinates(se_pix)

        # Modify the loops table and start res list for the system
        self.system.start_res_list[seid]=new_start_res
        self.system.add_se_to_loop_table(self.seid, new_start_res)
        self.system.update_system()

        return IMP.core.MonteCarloMoverResult([se.get_particle_index()], 1.0)

    def do_reject(self):
        se_pix = self.structure_elements[self.seid].get_particle_index()
       
        print("Rejecting")
        # Reset the loops table
        self.system.start_res_list[self.seid]=0
        self.system.remove_se_from_loop_table(self.seid)
        
        # Change the start_res for this SEID
        self.structure_elements[self.seid].set_start_res_key(0)
        
        # Zero out the coordinates we added
        self.zero_coordinates(se_pix)
        
        self.system.update_system()
        
class SSEDeletionMover(IMP.threading.SSEThreadMover):
    '''
    A move that proposes deletion of a random SE from the model
    '''
    def __init__(self, system, structure_elements):
        se_pis = [structure_elements[seid].get_particle_index() for seid in structure_elements.keys()]
        IMP.threading.SSEThreadMover.__init__(self, 
                system.model, 
                se_pis,
                system.sequence_hierarchy.get_particle_index())
        self.structure_elements = structure_elements
        self.system = system

    def get_se(self, seid):
        return self.structure_elements[seid]

    def do_propose(self):
        # Get a random structure element that is built
        rand_srs = list(self.structure_elements.keys())
        random.shuffle(rand_srs)

        for seid in rand_srs:
            se = self.structure_elements[seid]
            sr = se.get_start_res()
            if sr != 0:
                break

            # If none are zero, we can't add one.  Retnurn nothing.
            if seid == rand_srs[-1]:
                return IMP.core.MonteCarloMoverResult([se.get_particle_index()], 0.0)

        self.seid = seid
        self.old_start_res = se.get_start_res()
        se_pix = se.get_particle_index()

        # Change sequence coordinates to zero and set start_res_key of SE to zero
        self.zero_coordinates(se_pix)
        se.set_start_res_key(0)
        
        self.system.start_res_list[self.seid] = 0

        # Modify the loops table for the system
        self.system.remove_se_from_loop_table(self.seid)
        self.system.update_system()

        return IMP.core.MonteCarloMoverResult([se.get_particle_index()], 1.0)


    def do_reject(self):
        se = self.get_se(self.seid)
        se_pix = se.get_particle_index() 
        # Reset the loops table, reset start_Res_key and transform the coordinates back.
        self.system.add_se_to_loop_table(self.seid, self.old_start_res)
        se.set_start_res_key(self.old_start_res)
        self.transform_coordinates(se_pix)
        self.system.update_system()



class SSEShiftMover(IMP.threading.SSEThreadMover):
    '''
    A mover that shifts a random SSE start residue within the same loop
    '''
    def __init__(self, system, structure_elements, max_disp=None):
        se_pis = [structure_elements[seid].get_particle_index() for seid in structure_elements.keys()]
        #se_pis = [se.get_particle_index() for se in structure_elements]
        IMP.threading.SSEThreadMover.__init__(self,
                system.model,
                se_pis,
                system.sequence_hierarchy.get_particle_index())
        self.structure_elements = structure_elements
        self.system = system
        self.max_displacement = max_disp

    def get_se(self, seid):
        return self.structure_elements[seid]

    def do_propose(self):

        # Find random built SE
        rand_srs = list(self.structure_elements.keys())
        random.shuffle(rand_srs)
        
        for seid in rand_srs:
            se = self.structure_elements[seid]
            sr = se.get_start_res()
            if sr != 0:
                break

            # If none are built, we can't shift.  Retnurn nothing.
            if seid == rand_srs[-1]:
                self.seid=-1
                return IMP.core.MonteCarloMoverResult([se.get_particle_index()], 0.0)
        
        print("Shift", seid, sr)
        
        self.seid = seid
        self.old_start_res = sr
        se_pix = se.get_particle_index()

        # Zero the coordiantes
        self.zero_coordinates(se_pix)
        
        # find start residues within this loop
        start_res = self.system.get_available_start_residues(self.seid, exclude_self=True, same_loop=True)
        
        if len(start_res) == 0:
            print("ShiftMover WARNING: No available residues returned for this loop. Suspicious.")
            return IMP.core.MonteCarloMoverResult([se.get_particle_index()], 0.0)
        
        # Choose one of these at random
        new_start_res = random.choice(start_res)

        # Transform coordinates and update tables
        se.set_start_res_key(new_start_res)
        self.transform_coordinates(se_pix)
        self.system.start_res_list[self.seid]=new_start_res

        # Modify the loops table for the system
        self.system.modify_loop_table(self.seid, new_start_res)
        self.system.update_system()

        return IMP.core.MonteCarloMoverResult([se.get_particle_index()], 1.0)

    def do_reject(self):
        if self.seid == -1:
            return
        se = self.structure_elements[self.seid]
        se_pix = se.get_particle_index()
        new_sr = se.get_start_res()
        
        # Reset the loops table
        self.system.modify_loop_table(self.seid, self.old_start_res)
        self.system.start_res_list[self.seid]=self.old_start_res
        
        if new_sr != 0:
            self.zero_coordinates(se_pix)
        
        se.set_start_res_key(self.old_start_res)
        
        if self.old_start_res != 0:
            self.transform_coordinates(se_pix)
        self.system.update_system()

class SEMover(IMP.threading.StructureElementMover):
    # Mover wrapping an IMP.threading.StructureElementMover
    # This controls a single StructureElement
    def __init__(self, system, seid, chid = 0, zero_pct=50):
        IMP.threading.StructureElementMover.__init__(self, 
                system.model, 
                system.structure_elements[seid].get_particle_index(),
                system.sequence_hierarchy[chid].get_particle_index())
        print('particle indexes', system.structure_elements[seid].get_particle_index(), system.sequence_hierarchy[chid].get_particle_index())
        self.seid = seid
        self.system = system
        if zero_pct > 0:
            if zero_pct < 1.0:
                self.zero_pct = zero_pct*100
            elif zero_pct < 100:
                self.zero_pct = zero_pct
            else:
                raise Exception("zero_pct must be between zero and 100:", zero_pct)
        else:
            raise Exception("zero_pct must be between zero and 100:", zero_pct)

    def get_se(self):
        return self.system.structure_elements[self.seid]

    def do_propose(self):
        se = self.get_se()
        self.old_start_res = se.get_start_res()

        if self.old_start_res != 0:
            self.zero_coordinates()
        
        # if SE is not built, try to put it in:
        if self.old_start_res==0:
            start_res = self.system.get_available_start_residues(self.seid)
            if len(start_res)==0:
                return IMP.core.MonteCarloMoverResult([se.get_particle_index()], 0.0)

        # Otherwise, try to move it to another residue OR set it to zero 
        else:
            start_res = self.system.get_available_start_residues(self.seid, exclude_self=True)
            start_res = start_res+[0]*int(len(start_res)*self.zero_pct/100.0) # simple hack to add zeroes

        if len(start_res) == 0:
            return IMP.core.MonteCarloMoverResult([se.get_particle_index()], 0.0)
        new_start_res = random.choice(start_res)
        se.set_start_res_key(new_start_res)
        #print(self.system.get_built_residues())
        if new_start_res != 0:
            self.transform_coordinates()
        #print(self.system.get_built_residues())
        self.system.start_res_list[self.seid]=new_start_res

        # Modify the loops table for the system
        self.system.modify_loop_table(self.seid, new_start_res)

        return IMP.core.MonteCarloMoverResult([se.get_particle_index()], 1.0)



    def do_reject(self):
        new_sr = self.get_se().get_start_res()
        
        # Reset the loops table
        self.system.modify_loop_table(self.seid, self.old_start_res)
        
        if new_sr != 0:
            self.zero_coordinates()
        
        self.get_se().set_start_res_key(self.old_start_res)
        
        if self.old_start_res != 0:
            self.transform_coordinates()

def compare_threading_models(pdb_queries, pdb_reference, close_ca_threshold=5.0, exact_match_threshold=6.0):
    # Given two PDB files, return a number of metrics regarding the fit

    # Completeness:  How many residues built in query vs. reference

    # Close:  Number of query CA atoms within 5 angstroms of any reference CA atom. 

    # Sequence:  Number of query CA atoms within 5 angstroms of exact reference CA atom.

    if not hasattr(pdb_queries, "__iter__"):
        pdb_queries = pdb_queries

    m = IMP.Model()

    h = IMP.atom.read_pdb(pdb_reference, m, IMP.atom.CAlphaPDBSelector())

    ref_parts = IMP.atom.Selection(h, resolution=1).get_selected_particles()

    coords = {}
    for rp in ref_parts:
        resid = IMP.atom.Residue(IMP.atom.Hierarchy(rp).get_parent()).get_index()
        coords[resid] = IMP.core.XYZ(rp).get_coordinates()
        
    for pdb in pdb_queries:
        h2 = IMP.atom.read_pdb(pdb, m, IMP.atom.CAlphaPDBSelector())

        query_parts = IMP.atom.Selection(h2, resolution=1).get_selected_particles()
        
        query_coords = {}
        
        diffs = {}
        # Log coordinates and compute differences
        n_seq_close = 0
        n_ca_close = 0
        for pi in query_parts:
            resid = IMP.atom.Residue(IMP.atom.Hierarchy(pi).get_parent()).get_index()
            crds = IMP.core.XYZ(pi).get_coordinates()
            query_coords[resid] = crds
            try: 
                diffs[resid] = IMP.algebra.get_distance(crds, coords[resid]) 
                if diffs[resid] < exact_match_threshold:
                    n_seq_close+=1

                for ck in coords.keys():
                    if IMP.algebra.get_distance(crds, coords[ck]) < close_ca_threshold:
                        n_ca_close+=1
                        break

            except:
                continue

        print(len(coords.keys()), len(query_coords.keys()), n_seq_close, n_ca_close)


class SecondaryStructureRestraint(IMP.pmi.restraints.RestraintBase):

    """Restraint to keep all structures inside sphere."""

    def __init__(self,
                 hierarchy = None,
                 structural_elements=None,
                 prot = None, 
                 ss_file = None,
                 input_format = 'RaptorX',
                 weight = 1.0,
                 label = None):
        
        """Setup external barrier restraint.
        @param hierarchies IMP Hierarchy
        @param ss_file
        """

        
        self.ss_file = ss_file
        model = hierarchy.get_model()

        # Get particles in sequence chain
        
        
        super(SecondaryStructureRestraint, self).__init__(model,label=label,
                                                          weight=weight)

        if input_format == 'PSIPRED':
            self.read_psipred_ss2()
        elif input_format == 'RaptorX':
            self.read_raptorx_ss3()
            
        for ss in self.ss_probs:
            sp = IMP.atom.Selection(hierarchy,
                                    chain_id='S'+prot,
                                    residue_index = ss[0],
                                    copy_index=0).get_selected_particles()
            if len(sp)>0:
                IMP.atom.SecondaryStructureResidue.setup_particle(
                    sp[0], float(ss[1]), float(ss[2]), float(ss[3]))
                
        # Get SE particle indexes
        se_part_indexes = [s.get_particle_index() for s in structural_elements.values()]

        # Check number of particles
        #if len(ss_probs) != len(ps):
        #    raise Exception("SS file does not have same # of residues as protein")


        # Setup restraint
        chains = IMP.atom.get_by_type(hierarchy, IMP.atom.CHAIN_TYPE)
        seq_chain = [ch for ch in chains if 'seq_chain_'+prot in ch.get_name()][0]
        
        r = IMP.threading.SecondaryStructureParsimonyRestraint(model, se_part_indexes, seq_chain.get_particle_index(), 1.0)

        self.rs.add_restraint(r)


    def read_raptorx_ss3(self):
        '''
        IMP order: helix, strand, coil
        raptorX: helix, strand, coil
        '''
        flag = 0
        self.ss_probs = []
        for line in open(self.ss_file,'r'):
            if 'SS3' in line:
                flag = 1
            elif 'SS8' in line:
                flag = 0
                break
            if flag == 1:
                vals = line.split()
                if line[0] != '#' and len(vals)== 6:
                    self.ss_probs.append((int(vals[0]),float(vals[3]),float(vals[4]),float(vals[5])))
                    
    def read_psipred_ss2(self):
        '''
        IMP order: helix, strand, coil
        psi-pred order: helix, strand, coil
        '''

        f = open(self.ss_file, 'r')

        self.ss_probs = []
        for line in f.readlines():
            fields = line.strip().split()
            try: 
                int(fields[0])
                self.ss_probs.append((fields[4], fields[5], fields[3]))
            except:
                continue

    def get_name(self):
        return self.name

    def get_label(self):
        return self.label

class DistanceRestraint(IMP.pmi.restraints.RestraintBase):

    """A simple distance restraint"""

    def __init__(self,
                 representation=None,
                 tuple_selection1=None,
                 tuple_selection2=None,
                 distancemin=0,
                 distancemax=100,
                 resolution=1.0,
                 kappa=1.0,
                 root_hier=None,
                 label=None,
                 weight=1.):
        """Setup distance restraint.
        @param representation DEPRECATED
        @param tuple_selection1 (resnum, resnum, molecule name, copy
               number (=0))
        @param tuple_selection2 (resnum, resnum, molecule name, copy
               number (=0))
        @param distancemin The minimum dist
        @param distancemax The maximum dist
        @param resolution For selecting particles
        @param kappa The harmonic parameter
        @param root_hier The hierarchy to select from (use this instead of
               representation)
        @param label A unique label to be used in outputs and
                     particle/restraint names
        @param weight Weight of restraint
        @note Pass the same resnum twice to each tuple_selection. Optionally
              add a copy number (PMI2 only)
        """

        print('@@@@@@@@@@@@@@@@@@@@@@')
        
        if tuple_selection1 is None or tuple_selection2 is None:
            raise Exception("You must pass tuple_selection1/2")
        ts1 = IMP.core.HarmonicUpperBound(distancemax, kappa)
        ts2 = IMP.core.HarmonicLowerBound(distancemin, kappa)

        if representation and not root_hier:
            model = representation.prot.get_model()
            particles1 = IMP.pmi.tools.select(representation,
                                              resolution=resolution,
                                              name=tuple_selection1[2],
                                              residue=tuple_selection1[0])
            particles2 = IMP.pmi.tools.select(representation,
                                              resolution=resolution,
                                              name=tuple_selection2[2],
                                              residue=tuple_selection2[0])
        elif root_hier and not representation:
            model = root_hier.get_model()
            copy_num1 = 0
            if len(tuple_selection1) > 3:
                copy_num1 = tuple_selection1[3]
            copy_num2 = 0
            if len(tuple_selection2) > 3:
                copy_num2 = tuple_selection2[3]

            sel1 = IMP.atom.Selection(root_hier,
                                      resolution=resolution,
                                      chain_id=tuple_selection1[2],
                                      residue_index=tuple_selection1[0])
            particles1 = sel1.get_selected_particles()
            sel2 = IMP.atom.Selection(root_hier,
                                      resolution=resolution,
                                      chain_id=tuple_selection2[2],
                                      residue_index=tuple_selection2[0])
            particles2 = sel2.get_selected_particles()
        else:
            raise Exception("Pass representation or root_hier, not both")

        super(DistanceRestraint, self).__init__(model, label=label,
                                                weight=weight)

        print("Created distance restraint between "
              "%s and %s" % (particles1[0].get_name(),
                             particles2[0].get_name()))

        if len(particles1) > 1 or len(particles2) > 1:
            raise ValueError("more than one particle selected")

        self.rs.add_restraint(
            IMP.core.DistanceRestraint(self.model, ts1,
                                       particles1[0],
                                       particles2[0]))
        self.rs.add_restraint(
            IMP.core.DistanceRestraint(self.model, ts2,
                                       particles1[0],
                                       particles2[0]))

        
class ChainConnectivityRestraint_test(IMP.pmi.restraints.RestraintBase):

    def __init__(self,
                 root_hier,
                 chain,
                 slope = 1.0,
                 n_sds = 1.0,
                 label = None,
                 weight = 1.0):

        model = root_hier.get_model()
        super(ChainConnectivityRestraint_test, self).__init__(model, label=label, weight=weight)
        self.chain = chain
        self.n_sds = n_sds
        self.r = IMP.core.HarmonicUpperBound(0, slope)
        
        self.nhh_means = [3.81, 3.036, 2.836, 2.511, 2.275, 2.178, 2.026, 1.876, 1.835, 1.669, 1.658, 1.666, 1.625, 1.53,
                          1.445, 1.374, 1.292, 1.212, 1.164, 1.133, 1.049, 1.043, 1.074, 0.977, 0.965, 0.938, 0.868, 0.824, 0.805,
                          0.788]

        self.nhh_sds = [0.027, 0.284, 0.397, 0.441, 0.483, 0.499, 0.504, 0.537, 0.534, 0.538, 0.545, 0.507, 0.494, 0.468,
                        0.447, 0.428, 0.439, 0.415, 0.432, 0.392, 0.382, 0.38, 0.401, 0.381, 0.38, 0.317, 0.328, 0.304, 0.318,
                        0.273]

    def get_sd_distance_per_residue(self, nres):
        if nres > 29:
            sds = self.nhh_sds[29]
        else:
            sds = self.nhh_sds[nres]

        return sds

    def get_mean_distance_per_residue(self, nres):
        if nres > 29:
            mean = self.nhh_means[29]
        else:
            mean = self.nhh_means[nres]

        return mean
        
    def get_mean_distance(self, nres):
        return nres * self.get_mean_distance_per_residue(nres);

    def get_max_distance(self, nres):
        meandist = self.get_mean_distance(nres)
        sddist = self.get_sd_distance_per_residue(nres)
        #print('nres', nres, 'meandist', meandist, sddist, meandist + sddist * self.n_sds)
        return  meandist + sddist * self.n_sds;

    def evaluate(self):
        score = 0
        SEs_coors = []
        temp_XYZ = []
        for r in self.chain.get_children():
            resi_XYZ = IMP.core.XYZ(r)
            if resi_XYZ.get_coordinates_are_optimized():
                temp_XYZ.append([IMP.atom.Residue(r).get_index(), resi_XYZ])
            else:
                if len(temp_XYZ)>0:
                    SEs_coors.append(temp_XYZ)
                    temp_XYZ = []
                    
        if len(SEs_coors) == 0:
            return 0.0
        else:
            res_dif = [(SEs_coors[i+1][0][0]-SEs_coors[i][-1][0]) for i in range(len(SEs_coors)-1)]
            model_distances = [IMP.core.get_distance(SEs_coors[i+1][0][1],SEs_coors[i][-1][1]) for i in range(len(SEs_coors)-1)]
            max_distances = [self.get_max_distance(n) for n in res_dif]

            for d_mod, d_max, dd in zip(model_distances,max_distances, res_dif):
                d = d_mod-d_max
                #print('model distance', d_mod, 'd_max', d_max, dd)
                if d>0 and dd>0:
                    score += self.r.evaluate(d)
                elif dd< 0:
                    score += 1000
            return score
        

class ChainConnectivityRestraint(IMP.pmi.restraints.RestraintBase):

    def __init__(self,
                 root_hier,
                 SEs,
                 residues_N = None,
                 residues_C = None,
                 slope = 1.0,
                 n_sds = 2.0,
                 label = None,
                 weight = 1.0):

        model = root_hier.get_model()
        super(ChainConnectivityRestraint, self).__init__(model, label=label, weight=weight)
        self.SEs = SEs
        self.residues_N = residues_N
        self.residues_C = residues_C
        self.n_sds = n_sds
        
        self.r = IMP.core.HarmonicUpperBound(0, slope)
        
        self.nhh_means = [3.81, 3.036, 2.836, 2.511, 2.275, 2.178, 2.026, 1.876, 1.835, 1.669, 1.658, 1.666, 1.625, 1.53,
                          1.445, 1.374, 1.292, 1.212, 1.164, 1.133, 1.049, 1.043, 1.074, 0.977, 0.965, 0.938, 0.868, 0.824, 0.805,
                          0.788]

        self.nhh_sds = [0.027, 0.284, 0.397, 0.441, 0.483, 0.499, 0.504, 0.537, 0.534, 0.538, 0.545, 0.507, 0.494, 0.468,
                        0.447, 0.428, 0.439, 0.415, 0.432, 0.392, 0.382, 0.38, 0.401, 0.381, 0.38, 0.317, 0.328, 0.304, 0.318,
                        0.273]

    def get_sd_distance_per_residue(self, nres):
        if nres > 29:
            sds = self.nhh_sds[29]
        else:
            sds = self.nhh_sds[nres]

        return sds

    def get_mean_distance_per_residue(self, nres):
        if nres > 29:
            mean = self.nhh_means[29]
        else:
            mean = self.nhh_means[nres]

        return mean
        
    def get_mean_distance(self, nres):
        return nres * self.get_mean_distance_per_residue(nres);

    def get_max_distance(self, nres):
        meandist = self.get_mean_distance(nres)
        sddist = self.get_sd_distance_per_residue(nres)
        #print('nres', nres, 'meandist', meandist, sddist, self.n_sds, meandist + sddist * self.n_sds)
        return  meandist + sddist * self.n_sds

    def get_distance_fixed_end(self, chain_coords_sorted, chain_end_tag):

        if chain_end_tag == 'N':
            index = 0
            residues = self.residues_N
        elif chain_end_tag == 'C':
            index = -1
            residues = self.residues_C

        DD = []
        for p in residues:
            chain_p = IMP.atom.Hierarchy(p).get_parent().get_name()
            copy_p = IMP.atom.Copy(IMP.atom.Hierarchy(p).get_parent()).get_copy_index()
            name_p = f'{chain_p}_{copy_p}'
            if name_p in chain_coords_sorted.keys():
                SE_p = chain_coords_sorted[name_p][index][1]
                d_mod = IMP.algebra.get_distance(SE_p[index], IMP.core.XYZ(p).get_coordinates())
                n_res_dif = abs(IMP.atom.Residue(p).get_index() - chain_coords_sorted[name_p][index][0][index]) 
                d_max_t = self.get_max_distance(n_res_dif)
                d_max = self.lmax_flory(n_res_dif)
                d = d_mod - d_max
                #print('--------', d_mod, d_max, d) 
                DD.append(d)
        
        return DD

    def lmax_flory(self, n):
    	gamma = 6.046
    	delta = 3.46
    	if n %2. == 0:
        	return gamma * (n/2.-1) + delta
    	else:
        	return gamma * (n-1.)/2.    
 
    def evaluate(self):
        score = 0
        SEs_coors = []
        temp_XYZ = []

        # Get SE coordinates by chain
        chain_coords = {}
        for se in self.SEs:
            chain = se.get_chain()
            coords = se.get_all_coordinates()
            polarity = se.get_polarity()
            resi = numpy.arange(se.get_first_residue_number(),se.get_last_residue_number()+1,1)
            if polarity == -1:
                coords = coords[::-1]
            if chain in chain_coords.keys():
                chain_coords[chain].append([resi,coords])
            else:
                chain_coords[chain] = [[resi,coords]]
                
        # Sort the SE per chain
        chain_coords_sorted = {}
        for sel_chain, coors in chain_coords.items():
            temp_chain = sorted(chain_coords[sel_chain], key=lambda x:x[0][0])
            chain_coords_sorted[sel_chain] = temp_chain

        #print('chain_coords_sorted', chain_coords_sorted)
        # Evaluate restraint for terminal residues
        if self.residues_N:
            distances_N = self.get_distance_fixed_end(chain_coords_sorted, 'N')
            for d in distances_N:
                #print('N', d, self.r.evaluate(d))
                score += self.r.evaluate(d)

        if self.residues_C:
            distances_C = self.get_distance_fixed_end(chain_coords_sorted, 'C')
            for d in distances_C:
                #print('C', d, self.r.evaluate(d))
                score += self.r.evaluate(d)
        
        # Evaluate scores
        for sel_chain, temp_chain in chain_coords_sorted.items():
            res_dif = [temp_chain[i+1][0][0]-temp_chain[i][0][-1] for i in range(len(temp_chain)-1)]
            model_distances = [IMP.algebra.get_distance(temp_chain[i+1][1][0],temp_chain[i][1][-1]) for i in range(len(temp_chain)-1)]
            #print('model_distances', model_distances)
            max_distances = []
            for n in res_dif:
                if n<=12:
                   max_distances.append(self.get_max_distance(n))
                elif n>12 and n<=28:
                    d_s = self.get_max_distance(n)
                    d_f = self.lmax_flory(n)
                    w_f = (n-11)/(28-11)
                    w_s = 1-w_f 
                    max_distances.append(w_s*d_s + w_f*d_f)
                else:
                    max_distances.append(self.lmax_flory(n))
                
            max_distances_old = [self.get_max_distance(n) for n in res_dif]
            #print('distances SE', sel_chain, model_distances, max_distances, res_dif)
            #print(max_distances_old)
            
            for d_mod, d_max, dd in zip(model_distances,max_distances, res_dif):
                d = d_mod-d_max
                #print('model distance', d_mod, 'd_max', d_max, 'dd', dd)
                if d>0 and dd>0:
                    score += self.r.evaluate(d)
                    #print(self.r.evaluate(d))
                elif dd< 0:
                    #print('here!', 1000)
                    score += 10000
        return score
        

        
class StructuralElementConnectivityRestraint(IMP.pmi.restraints.RestraintBase):

    """A simple distance restraint"""

    def __init__(self,
                 root_hier,
                 SEs,
                 p1 = None,
                 p2 = None,
                 slope = 1.0,
                 dpr = 1.0,
                 label=None,
                 weight=1.):
        """Setup distance restraint.
        @param SEs list of all structural elements
        @param label A unique label to be used in outputs and
                     particle/restraint names
        @param weight Weight of restraint
        @note Pass the same resnum twice to each tuple_selection. Optionally
              add a copy number (PMI2 only)
        """

        model = root_hier.get_model()
        super(StructuralElementConnectivityRestraint, self).__init__(model, label=label, weight=weight)

        SEs_list  = list(SEs.values())
        for i, SE in enumerate(SEs_list[:-1]):
            p1 = SEs[i].get_particle_index()
            p2 = SEs[i+1].get_particle_index()
        
            r = IMP.threading.StructureElementConnectivityRestraint(model,
                                                                    IMP.core.HarmonicUpperBound(0, slope),
                                                                    p1,
                                                                    p2,
                                                                    dpr,
                                                                    "")

            #print(r.get_number_of_residues(), r.get_max_distance(), r.get_model_distance())
            
            self.rs.add_restraint(r)

    def get_label(self):
        return self.label

        '''
        for i, se in enumerate(SEs):
            p1 = SEs[i].get_particle_index()
            p2 = SEs[i+1].get_particle_index()

            r = IMP.threading.StructureElementConnectivityRestraint(model,
                                                                    IMP.core.HarmonicUpperBound(0, slope),
                                                                    p1,
                                                                    p2,
                                                                    dpr,
                                                                    "")
            self.rs.add_restraint(r)

            print(r.get_number_of_residues(), r.get_max_distance, r.get_model_distance())
            print("Created connectivity restraint between "
              "%s and %s" % (p1.get_name(),
                             p2.get_name()))

        '''

class StructuralElementCrossLinksRestraint(IMP.pmi.restraints.RestraintBase):

    """A simple distance restraint"""

    def __init__(self,
                 root_hier,
                 fname,
                 chain_mapping,
                 label=None,
                 weight=1.):
        
        """Setup distance restraint.
        @param SEs list of all structural elements
        @param label A unique label to be used in outputs and
                     particle/restraint names
        @param weight Weight of restraint
        @note Pass the same resnum twice to each tuple_selection. Optionally
              add a copy number (PMI2 only)
        """

        model = root_hier.get_model()
        super(StructuralElementCrossLinksRestraint, self).__init__(model, label=label, weight=weight)


        # Mapping protein to sequence chain
       
        missing = []
        nn = 0 
        for line in open(fname,'r'):
            vals = line.split(',')
            if vals[0]=='Protein1': continue
            prot1 = vals[0]
            prot2 = vals[2]
            r1 = vals[1]
            r2 = vals[3]
            
            xl_slope = 0.05
            xl_constant = 1.0
            length = 21.0

            mols = IMP.atom.get_by_type(root_hier, IMP.atom.MOLECULE_TYPE)
            mol_names = [m.get_name() for m in mols]
            if (prot1 =='Ig'):
                continue

            else:
                p1 = IMP.atom.Selection(root_hier, molecule=prot1, residue_index = int(r1)).get_selected_particles()
                p2 = IMP.atom.Selection(root_hier, molecule=prot2, residue_index = int(r2)).get_selected_particles()
                sp1 = IMP.atom.Selection(root_hier, chain_id='S'+prot1, residue_index = int(r1)).get_selected_particles()
                sp2 = IMP.atom.Selection(root_hier, chain_id='S'+prot2, residue_index = int(r2)).get_selected_particles()
                
                if len(p1)>0 and len(p2)>0:
                    #print(prot1, prot2, r1, r2, p1, p2, sp1, sp2)
                    r = IMP.threading.LoopPairDistanceRestraint(model, IMP.core.HarmonicUpperBound(length, xl_slope), 
                                                                p1, p2, xl_constant, length)
                    #self.rs.add_restraint(r)
                    #nn += 1
                elif len(sp1)>0 and len(p2)>0:
                    nn += 1
                    #print(prot1, prot2, r1, r2, p1, p2, sp1, sp2)
                    r = IMP.threading.LoopPairDistanceRestraint(model, IMP.core.HarmonicUpperBound(length, xl_slope), 
                                               sp1, p2, xl_constant, length)
                    self.rs.add_restraint(r)
                elif len(p1)>0 and len(sp2)>0:
                    nn += 1
                    #print(prot1, prot2, r1, r2, p1, p2, sp1, sp2)
                    r = IMP.threading.LoopPairDistanceRestraint(model, IMP.core.HarmonicUpperBound(length, xl_slope), 
                                               p1, sp2, xl_constant, length)
                    self.rs.add_restraint(r)
                elif len(sp1)>0 and len(sp2)>0:
                    nn += 1
                    #print(prot1, prot2, r1, r2, p1, p2, sp1, sp2)
                    r = IMP.threading.LoopPairDistanceRestraint(model, IMP.core.HarmonicUpperBound(length, xl_slope), 
                                               sp1, sp2, xl_constant, length)
                    self.rs.add_restraint(r)
        print('Effective number of XLs', nn)

    def get_label(self):
        return self.label
                

class SEDistanceRestraint(IMP.Restraint):

    """A simple distance restraint"""

    def __init__(self,
                 m,
                 func,
                 p1,p2,
                 name=None, weight=1.):
        
        #self.model = model

        self.p1 = p1
        self.p2 = p2
        self.func = func

        
        xl_constant = 1.0
        super(SEDistanceRestraint, self).__init__(m, name='test')

        length = 21.0
        xl_slope = 0.05
        self.xl_constant = 1.0
        self.A = IMP.threading.LoopPairDistanceRestraint(m,
                                                         IMP.core.HarmonicUpperBound(length, xl_slope),
                                                         p1, p2, self.xl_constant, length )
        

    def unprotected_evaluate(self, da):
        # Get coordinates
        
        XYZ1 = IMP.core.XYZ(self.p1)
        XYZ2 = IMP.core.XYZ(self.p2)

        is_coor1 = all(not np.isinf(c) and c>0 for c in XYZ1.get_coordinates())
        is_coor2 = all(not np.isinf(c) and c>0 for c in XYZ2.get_coordinates())
        score = 0.0
        if  is_coor1 and is_coor2:
            dist = IMP.core.get_distance(XYZ1, XYZ2)
            score = self.func.evaluate(dist)
        elif not is_coor1:
            # Find closets
            pn1 = self.A.get_closest_built_residue_particles(self.p1)
            if len(pn1) == 1:
                nres = abs(int(str(pn1[0].get_index()))-int(str(self.p1.get_index())))
                XYZn1 = IMP.core.XYZ(pn1[0])
                da = IMP.core.get_distance(XYZn1, XYZ2)
                db = nres * (nhh_mean[nres] + hh_std[nres]*self.xl_constant)
                dist = da + db
                score = self.func.evaluate(dist)
            elif len(pn1) == 2:
                print('Missing')
        elif not is_coor1:
            pn2 = self.A.get_closest_built_residue_particles(self.p2)
            if len(pn2) == 1:
                nres = abs(int(str(pn2[0].get_index()))-int(str(self.p2.get_index())))
                XYZn2 = IMP.core.XYZ(pn2[0])
                da = IMP.core.get_distance(XYZn2, XYZ1)
                db = nres * (nhh_mean[nres] + hh_std[nres]*self.xl_constant)
                dist = da + db
                score = self.func.evaluate(dist)
            elif len(pn2) == 2:
                print('Missing')
                
        return score 
    

            
