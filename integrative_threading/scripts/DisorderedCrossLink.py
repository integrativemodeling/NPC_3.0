import numpy
import math
import ihm.restraint


nsh_mean = [0, 3.809, 3.137, 2.818, 2.482, 2.154, 1.928, 1.749, 1.67, 1.531, 1.428, 1.377, 1.282, 1.261, 1.203,
                1.135, 1.045, 1.004, 1.02, 0.977, 0.928, 0.865, 0.834, 0.811, 0.756, 0.761, 0.749, 0.777, 0.74, 0.655,
                0.648]
nhs_mean = [0, 3.809, 3.137, 2.818, 2.482, 2.154, 1.928, 1.749, 1.67, 1.531, 1.428, 1.377, 1.282, 1.261, 1.203,
                1.135, 1.045, 1.004, 1.02, 0.977, 0.928, 0.865, 0.834, 0.811, 0.756, 0.761, 0.749, 0.777, 0.74, 0.655,
                0.648]
nhh_mean = [0, 3.81, 3.036, 2.836, 2.511, 2.275, 2.178, 2.026, 1.876, 1.835, 1.669, 1.658, 1.666, 1.625, 1.53,
                1.445, 1.374, 1.292, 1.212, 1.164, 1.133, 1.049, 1.043, 1.074, 0.977, 0.965, 0.938, 0.868, 0.824, 0.805,
                0.788]
nss_mean = [0, 3.81, 3.19, 1.846, 1.607, 1.274, 1.14, 1.139, 1.198, 1.177, 1.115, 1.029, 1.048, 0.935, 0.91, 0.908,
                0.85, 0.83, 0.852, 0.849, 0.761, 0.722, 0.742, 0.684, 0.677, 0.611, 0.587, 0.596, 0.565, 0.576, 0.532]


# Standard deviations of the above normalised distances

hh_std = [0, 0.027, 0.284, 0.397, 0.441, 0.483, 0.499, 0.504, 0.537, 0.534, 0.538, 0.545, 0.507, 0.494, 0.468,
              0.447, 0.428, 0.439, 0.415, 0.432, 0.392, 0.382, 0.38, 0.401, 0.381, 0.38, 0.317, 0.328, 0.304, 0.318,
              0.273]
ss_std = [0, 0.027, 0.313, 0.293, 0.469, 0.419, 0.474, 0.49, 0.505, 0.447, 0.501, 0.475, 0.479, 0.417, 0.451, 0.416,
              0.373, 0.395, 0.47, 0.418, 0.36, 0.349, 0.359, 0.312, 0.302, 0.281, 0.279, 0.264, 0.259, 0.346, 0.257]
sh_std = [0, 0.067, 0.278, 0.361, 0.418, 0.45, 0.448, 0.455, 0.436, 0.452, 0.438, 0.416, 0.407, 0.402, 0.411, 0.405,
              0.381, 0.378, 0.373, 0.36, 0.372, 0.338, 0.322, 0.308, 0.285, 0.289, 0.296, 0.298, 0.294, 0.286, 0.208]
hs_std = [0, 0.067, 0.278, 0.361, 0.418, 0.45, 0.448, 0.455, 0.436, 0.452, 0.438, 0.416, 0.407, 0.402, 0.411, 0.405,
              0.381, 0.378, 0.373, 0.36, 0.372, 0.338, 0.322, 0.308, 0.285, 0.289, 0.296, 0.298, 0.294, 0.286, 0.208]


pi_over_2 = math.pi * 0.5

def get_closest_endpoints(res, atoms):
    # returns two tuples of: ((x, y, z), nres), one for the c-terminal, one
    # for the n-terminal side.
    #

    res_dict = {}

    for a in atoms:
        res_dict[a.seq_id_range[0]] = a

    if res in res_dict.keys():
        a = res_dict[res]
        xyz = numpy.array([a.x, a.y, a.z])
        #print(res, res)
        #print((xyz, 0), (xyz, 0))
        return (xyz, 0), (xyz, 0)

    else:
        for i in range(1,res):
            if res-i in res_dict.keys():
                if i > 30:
                    j=30
                else:
                    j=i
                a = res_dict[res-i]
                n_xyz=(numpy.array([a.x, a.y, a.z]), i*nhs_mean[j])
                #print(res, res-i)
                break

            elif i==res-1:
                if i > 30:
                    j=30
                else:
                    j=i
                n_xyz = (numpy.array([0,0,0]), i*nhs_mean[j])
                #print(res, res-i)
                break

        for i in range(1,res):
            
            if res+i in res_dict.keys():

                if i > 30:
                    j=30
                else:
                    j=i
                a = res_dict[res+i]
                c_xyz=(numpy.array([a.x, a.y, a.z]), i*nhs_mean[j])
                #print(res, res+i)
                break

            elif i==res+1:
                if i > 30:
                    j=30
                else:
                    j=i
                c_xyz = (numpy.array([0,0,0]), i*nhs_mean[j])
                #print(res, res+i)
                break
    #print("NXYZ", n_xyz, c_xyz)
    #print("D_vec", n_xyz[0]-c_xyz[0])
    return n_xyz, c_xyz

def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))

def length(v):
  return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
    #print("DOT", dotproduct(v1, v2),length(v1), length(v2), dotproduct(v1, v2) / (length(v1) * length(v2)))
    ip = dotproduct(v1, v2) / (length(v1) * length(v2))
    if ip < 1.0:
        return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))
    else:
        return 0.0

def calculate_sphere_cap_distance(R, d, alpha):
    sina = math.sin(alpha);
    #print("CALC", R, d, alpha, sina, d*d*(sina*sina - 1) + R*R)
    dist = -d * sina + math.sqrt(d*d*(sina*sina - 1) + R*R)
    if (dist < 0): 
        dist = -1*dist
    return dist


def calculate_pseudosites(n1, c1, n2, c2):

    # First, get the centers and distances
    pdist1 = length(n1[0]-c1[0])
    pdist2 = length(n2[0]-c2[0])

    pcen1 = (n1[0]+c1[0])/2
    pcen2 = (n2[0]+c2[0])/2

    #print(" -- P", pcen1-pcen2, dotproduct(pcen1-pcen2, n2[0]-c2[0]))

    xl_vec = pcen1-pcen2
    xl_dist = length(xl_vec)

    if xl_dist == 0:
        return pcen1, pcen2

    offset_1 = 0
    offset_2 = 0

    if pdist1 > 0:
        #print("pdist1", pdist1)
        D_vec_1 = n1[0]-c1[0]

        angle1 = angle(D_vec_1, xl_vec)

        if angle1 > pi_over_2:
            angle1 = angle1 - pi_over_2
            cen_dist1 = length(pcen1-c1[0])
        else:
            cen_dist1 = length(pcen1-n1[0])

        # This is the case where the spheres do not intersect
        if length(D_vec_1) > n1[1] + c1[1]:
            offset_1 = 0
        elif n1[1] > c1[1] + length(D_vec_1):
            offset_1 = n1[1]
        elif c1[1] > n1[1] + length(D_vec_1):
            offset_1 = c1[1]
        else:
            offset_1 = calculate_sphere_cap_distance(c1[1], cen_dist1, angle1);

    if pdist2 > 0:
        #print("pdist2", pdist2, n2[0]-c2[0], xl_vec)
        #print("                  ", pcen1, pcen2)
        D_vec_2 = n2[0]-c2[0]
        angle2 = angle(D_vec_2, xl_vec)

        if angle2 > pi_over_2:
            angle2 = angle2 - pi_over_2
            cen_dist2 = length(pcen2-n2[0])
        else:
            cen_dist2 = length(pcen2-c2[0])
        # This is the case where the spheres do not intersect
        if length(D_vec_2) > n2[1] + c2[1]:
            offset_2 = 0
        elif n1[1] > c2[1] + length(D_vec_2):
            offset_2 = n2[1]
        elif c1[1] > n2[1] + length(D_vec_2):
            offset_2 = c2[1]
        else:
            offset_2 = calculate_sphere_cap_distance(c2[1], cen_dist2, angle2);

    unit_vector = numpy.linalg.norm(xl_vec)

    # Psuedosites are along xl_vec with magnitude = offset
    ps1 = pcen1 - unit_vector * offset_1
    ps2 = pcen2 + unit_vector * offset_2

    #print("  -", offset_1, " | ", offset_2)

    return ps1, ps2


class PseudoSiteFinder(object):
    def __init__(self, model):
        self.model = model
        self.spheres = model._spheres
        # dict with keys=residues, values=PseudoSiteFeature
        self._seen_sites = {}

    def get_pseudo_site_endpoints(self, res1, res2):
        #print "*****************************"
        nte1, cte1 = get_closest_endpoints(res1, self.spheres)
        nte2, cte2 = get_closest_endpoints(res2, self.spheres)

        # If residue number is zero, this is a structured residue
        # The pseudo-sites are the coordinates of this residue
        ps1, ps2 = calculate_pseudosites(nte1, cte1, nte2, cte2)

        form = 0
        if nte1[1] == 0:
            form+=1
        if nte2[1] == 0:
            form+=2

        return (self.get_pseudo_site(ps1, res1),
                self.get_pseudo_site(ps2, res2), form)

    def get_pseudo_site(self, coord, res):
        if res not in self._seen_sites:
            # Don't duplicate features if we see the same residue multiple times
            ps = ihm.restraint.PseudoSite(
                    x=coord[0], y=coord[1], z=coord[2],
                    description='This pseudo site corresponds to residue '
                        + str(res) + ' from entity 1. It is not part of the '
                        'model representation but is part of the '
                        'experimental restraints')
            psf = ihm.restraint.CrossLinkPseudoSite(ps, model=self.model)
            self._seen_sites[res] = psf
        return self._seen_sites[res]
