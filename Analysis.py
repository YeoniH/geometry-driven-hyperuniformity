# Author: Sungyeon Hong (unless stated otherwise; other contributors are annotated)
# The "Analysis" module provides functions for a large set of computational analyses
# of geometry-driven hyperuniformity obtained through Lloyd relaxation/iteration process. 
# The cleaning of this lengthy module is in progress, so please feel free to contact 
# the author when confusion arises. 
# Also, you can find more details in "README.md".

from os import system, path
import time
import math
import numpy as np
import matplotlib.pyplot as plt
import freud
import collections
import itertools
from scipy.spatial import ConvexHull
from collections import defaultdict


class ConvexGeometry(object):
    
    def __init__(self, vertices):
        self.vertices = vertices

    def area(self):
        'Area of cross-section.'
        vertices = self.vertices
        vertices = list(vertices)
        # if vertices[0] != vertices[-1]:
        vertices = vertices + vertices[:1]
        x = [c[0] for c in vertices]
        y = [c[1] for c in vertices]
        s = 0
        for i in range(len(vertices) - 1):
            s += x[i] * y[i + 1] - x[i + 1] * y[i]
        A = s / 2
        return A
        
    def centroid(self):
        'Location of centroid.'
        # if vertices[0] != vertices[-1]:
        #     vertices = vertices + vertices[:1]
        vertices = self.vertices
        vertices = list(vertices)
        vertices = vertices + vertices[:1]
        x = [c[0] for c in vertices]
        y = [c[1] for c in vertices]
        sx = sy = 0
        A = self.area()
        for i in range(len(vertices) - 1):
            sx += (x[i] + x[i + 1]) * (x[i] * y[i + 1] - x[i + 1] * y[i])
            sy += (y[i] + y[i + 1]) * (x[i] * y[i + 1] - x[i + 1] * y[i])
        Cx = sx / (6 * A)
        Cy = sy / (6 * A)
        return np.array([float(f"{Cx:.16f}"), float(f"{Cy:.16f}"), 0.0])

    def quantizer_error(self, generating_point):
        'Moments and product of inertia about centroid.'
        vertices = self.vertices
        vertices = list(vertices)
        vertices = vertices + vertices[:1]
        A = self.area()
        Px, Py = generating_point[:2]
        Cx, Cy = self.centroid()[:2]
        x = [c[0] - Px for c in vertices]  # centred at P (generating point)
        y = [c[1] - Py for c in vertices]
        sxx = syy = sxy = 0
        for i in range(len(vertices) - 1):
            sxx += (y[i] ** 2 + y[i] * y[i + 1] + y[i + 1] ** 2) * (x[i] * y[i + 1] - x[i + 1] * y[i])
            syy += (x[i] ** 2 + x[i] * x[i + 1] + x[i + 1] ** 2) * (x[i] * y[i + 1] - x[i + 1] * y[i])
            sxy += (x[i] * y[i + 1] + 2 * x[i] * y[i] + 2 * x[i + 1] * y[i + 1] + x[i + 1] * y[i]) \
                   * (x[i] * y[i + 1] - x[i + 1] * y[i])
        Ixx = sxx / 12  # - A * Py ** 2    # Parallel axis theorem used when centred at P (generating point)
        Iyy = syy / 12  # - A * Px ** 2
        Ixy = sxy / 24  # - A * Px * Py
        return Ixx, Iyy, Ixy
    
class PeriodicVoro(object):
    def __init__(self, dim, points, L, run, step):
        self.dim = dim
        self.points = points
        self.num_points = len(points)
        self.L = L
        self.run = run
        self.step = step

        self.box = freud.box.Box.square(L)
        voro = freud.locality.Voronoi()
        voro.compute((self.box, points))
        self.voro = voro

    def update_by_lloyd_algorithm(self, pbc=True):
        L = self.L
        updated_points = np.zeros(self.points.shape)
        for i, polygon in enumerate(self.voro.polytopes):
            centroid = ConvexGeometry(polygon).centroid()
            x, y = centroid[0], centroid[1]
            if pbc:
                if x < -L / 2:
                    x += L
                elif x > L / 2:
                    x -= L
                if y < -L / 2:
                    y += L
                elif y > L / 2:
                    y -= L
            updated_points[i] = [x, y, 0.]
        return np.array(updated_points)

class DynamicalAnalysis(object):
    
    def __init__(self, coord, name, step, data_path, save_path, boxsz=None):

        self.coord = coord
        self.numPoints = len(coord)
        self.name = name
        self.step = step
        self.fullName = '{}-step-{}'.format(name, step)
        self.dataPath = data_path
        self.savePath = save_path
        
        fig = plt.figure(figsize=(5, 3))
        self.ax = fig.add_subplot(111)
        
        if boxsz == None:
            self.L_x = np.sqrt(self.numPoints)
            self.L_y = np.sqrt(self.numPoints)
        elif len(boxsz) == 1:
            self.L_x = boxsz
            self.L_y = boxsz
        elif len(boxsz) == 2:
            self.L_x = boxsz[0]
            self.L_y = boxsz[1]
        
        self.box = freud.box.Box.from_box([self.L_x, self.L_y])
        
        self.points = np.c_[self.coord, [0] * self.numPoints]
        self.voro = freud.locality.Voronoi()
        self.voro.compute((self.box, self.points))
        self.cluster = freud.cluster.Cluster()
        self.cluster.compute((self.box, self.points), )

        self.points_next = np.c_[coord_next, [0] * self.numPoints]
        self.voro_next = freud.locality.Voronoi()
        self.voro_next.compute((self.box, self.points_next))
    
    def msd(self):
        """ This function computes the meas-square displacement of a system """
        msd = freud.msd.MSD(box=self.box, mode="direct")
        msd.compute(positions=(self.points, self.points_next))
        msd.plot(ax=self.ax)

    def quantizer_energy(self):
        energy = []
        for idx, polytope in enumerate(self.voro.polytopes):
            vertices = polytope[:, :2]
            inertia = ConvexGeometry(vertices).quantizer_error(self.points[idx])
            quantizer_energy = np.abs(sum(inertia[:2]))
            energy.append(quantizer_energy)
        return energy
        
    def compute_total_quantizer_energy(self):
        quantizer_energy = self.quantizer_energy_of_voronoi_cells()
        total_quantizer_energy = float((self.numPoints**(2/self.dim)) / (self.dim*(self.L_x*self.L_y)**(1+2/self.dim))) * float(sum(quantizer_energy))
        return total_quantizer_energy


    ###### HEXATIC ORDER PARAMETER
    
    def compute_hexatic_order(self):
        op = freud.order.Hexatic(k=6)
        op.compute(system=({'Lx': self.L, 'Ly': self.L, 'dimensions': 2}, self.points))
        return op.particle_order

    def compute_hexatic_order_modulus(self):
        op = freud.order.Hexatic(k=6)
        op.compute(system=({'Lx': self.L, 'Ly': self.L, 'dimensions': 2}, self.points))
        return np.absolute(op.particle_order)

    def compute_hexatic_order_phase_angle(self, deg=False):
        op = freud.order.Hexatic(k=6)
        op.compute(system=({'Lx': self.L, 'Ly': self.L, 'dimensions': 2}, self.points))
        if deg:
            return np.angle(op.particle_order, deg=True)
        else:
            return np.angle(op.particle_order)  # in radians
            
    def compute_translational_order_modulus(self):
        """
        The translational order parameter in a two-dimensional crystal refers to the degree of
        positional correlation between particles in the crystal lattice.
        The translational order parameter is determined by analyzing the decay behavior of the
        translational correlation function, which reflects the symmetry and arrangement of
        particles in the crystal lattice.
        """
        top = freud.order.Translational(k=6)
        top.compute(system=({'Lx': self.L, 'Ly': self.L, 'dimensions': 2}, self.points))
        return np.absolute(top.particle_order)


    ###### CORRELATION FUNCTION
    
    def compute_correlation_function(self, to_correlate="phase_angle", bins=25):
        if to_correlate == "hex_order":
            order_param = self.compute_hexatic_order_modulus()
            # ylabel = r'$c_{|\psi_{6}|}(r)$'
        elif to_correlate == "bond_order":
            order_param = self.compute_hexatic_order()
            # ylabel = r'$c_{\psi_{6}}(r)$'
        elif to_correlate == "phase_angle":
            order_param = self.compute_hexatic_order_phase_angle()
            # ylabel = r'$c_{\theta}(r)$'

        mean_order_param = np.mean(order_param)
        values = []  # phase angle fluctuations
        for i in range(self.num_points):
            values.append(order_param[i] - mean_order_param)

        """ Create a CorrelationFunction object """
        cf = freud.density.CorrelationFunction(bins=bins, r_max=self.box.Lx / 2.01)
        cf.compute(system=(self.box, self.points), values=values)
        autocorr = np.c_[cf.bin_centers, cf.correlation]  # np.c_[np.abs(cf.bin_centers), np.abs(cf.correlation)]
        return autocorr

    def correlation_length(self, to_correlate="phase_angle", bins=25, norm=True, from_bin=1):
        autocorr = self.compute_correlation_function(to_correlate=to_correlate, bins=bins)[from_bin:]
        x_dist = autocorr.T[0]
        if norm:
            y_corr = autocorr.T[1]/autocorr.T[1][0]     # normalised!
        else:
            y_corr = autocorr.T[1]
        c_zero_at = 0
        for i in range(1, len(y_corr)):
            if y_corr[i - 1]*y_corr[i] < 0:
                a = (y_corr[i] - y_corr[i - 1]) / (x_dist[i] - x_dist[i - 1])
                c_zero_at = (1 / a) * (a * x_dist[i] - y_corr[i])
                break
        corr_length = c_zero_at
        return corr_length

    
    ###### T1 TRANSITION
    
    def group_quadruples(self):
        nb_dict = defaultdict(list)
        for n_i, nb in enumerate(self.voro.nlist):
            ii, jj = nb
            nb_dict[ii].append(jj)
        # print(nb_dict)

        quadruples = []
        sets_of_quadruples = []
        for i in range(self.numPoints):
            neighbors = nb_dict[i]
            pairs = list(itertools.combinations(neighbors, 2))
            for pair in pairs:
                quadruple_i = [i]
                pair = np.array(pair)
                if pair[1] in nb_dict[pair[0]]:
                    quadruple_i.append(pair[0])
                    quadruple_i.append(pair[1])
                if len(quadruple_i) == 3:
                    common_nb = list(set(nb_dict[quadruple_i[-1]]).intersection(nb_dict[quadruple_i[-2]]))
                    for idx in common_nb:
                        if idx != i:
                            quadruple_i.append(idx)
                    if set(quadruple_i) not in sets_of_quadruples:
                        sets_of_quadruples.append(set(quadruple_i))
                        quadruples.append(quadruple_i)
        return quadruples

    def t1_active_quadruples(self, quadruples, quadruples_next, empty_save=False):
        print("Extracting T1-quadruples:", self.fullName)
        start_time = time.time()
        Q_t1active = []
        edge_ratio_t1 = []
        for quad in quadruples:
            if len(quad) == 4:
                q0, q1, q2, q3 = quad
                old_link = []
                new_link = []
                if [q1, q0, q3, q2] in quadruples_next:
                    old_link = [q1, q2]
                    new_link = [q0, q3]
                elif [q2, q0, q3, q1] in quadruples_next:
                    old_link = [q1, q2]
                    new_link = [q0, q3]
                elif [q1, q3, q0, q2] in quadruples_next:
                    old_link = [q1, q2]
                    new_link = [q0, q3]
                elif [q2, q3, q0, q1] in quadruples_next:
                    old_link = [q1, q2]
                    new_link = [q0, q3]
                if old_link + new_link:
                    # t1_quad = old_link + new_link
                    Q_t1active.append(old_link + new_link)
            else:
                print(self.fullName, ', T1-active cells for non-quadruples:', quad)
        flatten_t1Q = [elt for quad in Q_t1active for elt in quad]
        if empty_save:   # Save empty files
            np.savetxt('{0}/t1_transition/T1-quad/{1}/T1-quad_{2}_t{3}.txt'.format(self.savePath, self.name, self.name[3:], self.step),
                       np.reshape(Q_t1active, (len(Q_t1active), 4)), fmt='%i %i %i %i')
        else:           # Don't save empty files
            if Q_t1active:
                np.savetxt('{0}/t1_transition/T1-quad/{1}/T1-quad_{2}_t{3}.txt'.format(self.savePath, self.name, self.name[3:],self.step),
                           np.reshape(Q_t1active, (len(Q_t1active), 4)), fmt='%i %i %i %i')
        print('Time taken:', time.strftime('%H:%M:%S', time.gmtime(time.time() - start_time)))
        return flatten_t1Q

    def compute_polygonal_symmetry(self):
        k_dict = {}
        for idx, polytope in enumerate(self.voro.polytopes):
            k_dict[idx] = len(polytope)
        magnitude = []
        angle = []
        for i, k in enumerate(k_dict.values()):
            op = freud.order.Hexatic(k=k)
            op.compute(system=({'Lx': self.L_x, 'Ly': self.L_y, 'dimensions': 2}, self.points),
                       neighbors={'num_neighbors': k})
            magnitude.append(np.absolute(op.particle_order[i]))
            angle.append(np.angle(op.particle_order[i]))
        return magnitude, angle

    def mean_polygonal_symmetry_of_k_skeleton(self, k=2):
        f_skeletons = open(
            '{0}/defects/{1}/skeletons/skeletons_{2}.txt'.format(self.savePath, self.name, self.fullName), 'r')
        k_skeletons = []
        for line in f_skeletons.readlines():
            line = line.strip().split(' ')
            if len(line) == k:
                k_skeletons.append([int(line[j]) for j in range(len(line))])
        symm_order_all = self.compute_polygonal_symmetry(order="magnitude")
        symm_order_k_skeletons = []
        for skeleton in k_skeletons:
            for elt in skeleton:
                symm_order_k_skeletons.append(symm_order_all[elt])
        return np.mean(symm_order_k_skeletons)

    
    ###### ISOPERIMETRIC RATIO

    def isoperimetric_ratios(self):
        ratio = []
        for idx, polytope in enumerate(self.voro.polytopes):
            vertices = polytope[:, :2]
            hull = ConvexHull(vertices)  # [1:])
            P = hull.area  # surface area
            A = hull.volume  # bounded volume
            isoperimetric_ratio = 4 * np.pi * float(A / P ** 2)
            ratio.append(isoperimetric_ratio)
        return ratio


    ###### TOPOLOGICAL CHARGE

    def topological_charge_of_every_polygon(self):
        charges = []
        for idx, polytope in enumerate(self.voro.polytopes):
            charge.append(int(6 - len(polytope)))
        return charges


    ###### DEFECT ANALYSIS
    
    def extract_defect_skeletons(self):
        defective_cells = {}
        for idx, polytope in enumerate(self.voro.polytopes):
            if len(polytope) in [5, 7]:
                defective_cells[idx] = len(polytope)
        # print(defective_cells.keys())

        nlist_defects = defaultdict(list)
        count = 0
        for idx, polytope in enumerate(self.voro.polytopes):
            # count += len(polytope)
            if idx not in defective_cells:
                count += len(polytope)
            else:
                assert len(polytope) == defective_cells[idx]
                for idx_neighbor in range(len(polytope)):
                    ii, jj = self.voro.nlist[idx_neighbor + count]
                    if ii not in defective_cells.keys():
                        print(ii)
                        # nlist_defects[ii].append(None)
                        count += 1
                    else:
                        if jj in defective_cells.keys():
                            nlist_defects[ii].append(jj)
                # print(idx, len(polytope))
                # print(self.voro.nlist[count:count+len(polytope)])
                count += len(polytope)
        isolated_defects = sorted(defective_cells.keys() - nlist_defects.keys())
        
        system("mkdir {0}/defects/{1} {0}/defects/{1}/nlist {0}/defects/{1}/skeletons {0}/defects/{1}/charge".format(self.savePath, self.name))
        f_nlist = open('{0}/defects/{1}/nlist/nlist_{2}.txt'.format(self.savePath, self.name, self.fullName), 'w')
        for defect in defective_cells.keys():
            if defect not in isolated_defects:
                for neighbor in nlist_defects[defect]:
                    f_nlist.write(str(defect)+' '+str(neighbor)+'\n')

        endpoints = []
        for defect in defective_cells.keys():
            if len(nlist_defects[defect]) == 1:
                endpoints.append(defect)
                
        print('Extracting defect skeletons with topo charge for ', self.fullName)
        already_grouped = []
        another_end = []
        clusters = []
        for p0 in endpoints:
            if p0 not in already_grouped:# and p0 not in another_end:
                # print('Start:', p0)#, nlist_defects[p0])
                chain = [p0]
                n1 = nlist_defects[p0][0]
                chain.append(n1)
                already_grouped = list(set(already_grouped + chain))
                while chain[-1] not in endpoints:
                    previous_chain_size = len(chain)
                    for n2 in nlist_defects[chain[-1]]:
                        if n2 not in chain:
                            chain.append(n2)
                    for item in chain:
                        for nn in nlist_defects[item]:
                            if nn not in chain:
                                chain.append(nn)
                    already_grouped = list(set(already_grouped + chain))
                    if len(chain) - previous_chain_size == 0:  # No update!
                        break
                if len(chain) != 0:
                    another_end.append(chain[-1])

                clusters.append(chain)

        for item in isolated_defects:
            clusters.append([item])

        already_sorted = []
        for item in defective_cells.keys():
            if item not in already_sorted and item not in already_grouped + isolated_defects:
                ring = [item]
                for n1 in nlist_defects[item]:
                    ring.append(n1)
                for member in ring:
                    for nn in nlist_defects[member]:
                        if nn not in ring:
                            ring.append(nn)
                already_sorted = list(set(already_sorted + ring))
                clusters.append(ring)

        chain_size = []
        for labeled_chain in clusters:
            chain_size.append(len(labeled_chain))
        print('Max chain size of {}: {}'.format(self.fullName, np.max(chain_size)))

        f_skeletons = open('{0}/defects/{1}/skeletons/skeletons_{2}.txt'.format(self.savePath, self.name, self.fullName), 'w')
        f_charge = open('{0}/defects/{1}/charge/charge_{2}.txt'.format(self.savePath, self.name, self.fullName), 'w')
        k_list = []
        for independent_chain in clusters:
            topological_charge = 0
            skeleton = []
            for item in independent_chain:
                f_skeletons.write(str(item) + ' ')
                skeleton.append(int(item))
                if defective_cells[item] == 5:
                    topological_charge += 1
                elif defective_cells[item] == 7:
                    topological_charge -= 1
                else:
                    topological_charge += 0
            f_charge.write(str(topological_charge) + ' \n')
            k_list.append(skeleton)
            if independent_chain:
                f_skeletons.write('\n')
        f_skeletons.close()
        f_charge.close()

        # f_skeletons = open('{0}/defects/{1}/skeletons/skeletons_{2}.txt'.format(self.savePath, self.name, self.fullName), 'r')
        # k_list = []
        # for i, line in enumerate(f_skeletons.readlines()):
        #     line = line.strip().split(' ')
        #     k_list.append(len(line))
        # k_list = list(set(k_list))
        # f_skeletons.close()

        size_sum = 0
        all_in_chains = []
        for independent_chain in clusters:
            size_sum += len(independent_chain)
            for item in independent_chain:
                all_in_chains.append(item)
        print('Total number of (labelled) defects of {}: {} ({})'.format(self.fullName, size_sum,
                                                                         len(defective_cells)))
        if size_sum != len(defective_cells):
            print('There must be duplicated cells in {}!'.format(self.fullName))

        return k_list

    def extract_skins_of_defect_skeletons(self):
        defective_cells = {}
        nlist_defect_skin = defaultdict(list)
        count = 0
        for idx, polytope in enumerate(self.voro.polytopes):
            if len(polytope) not in [5, 7]:
                count += len(polytope)
            else:
                defective_cells[idx] = len(polytope)
                for idx_neighbor in range(len(polytope)):
                    ii, jj = self.voro.nlist[idx_neighbor + count]
                    if ii not in defective_cells.keys():
                        count += 1
                    else:
                        # if jj not in defective_cells.keys():
                        if len(self.voro.polytopes[jj]) == 6:
                            nlist_defect_skin[ii].append(jj)
                            # print(len(self.voro.polytopes[jj]))
                count += len(polytope)

        system('mkdir {0}/defects/{1}/bodies'.format(self.savePath, self.name))
        f_skeletons = open('{0}/defects/{1}/skeletons/skeletons_{2}.txt'.format(self.savePath, self.name, self.fullName), 'r')
        f_skins = open('{0}/defects/{1}/bodies/skins_{2}.txt'.format(self.savePath, self.name, self.fullName), 'w')
        skins = []
        for i, line in enumerate(f_skeletons.readlines()):
            line = line.strip().split(' ')
            line = [int(line[j]) for j in range(len(line))]
            skins = []
            added_already = []
            for defect in line:
                skins.append(nlist_defect_skin[defect])
                for first_nb in nlist_defect_skin[defect]:
                    if first_nb not in added_already:
                        f_skins.write('{} '.format(first_nb))
                        added_already.append(first_nb)
            f_skins.write('\n')
            skins.append(skins)
        # print(skins)
        f_skins.close()
        print('Done with extracting defect skins in {} !'.format(self.fullName))

    def extract_independent_defect_bodies(self):
        print('Extracting defect bodies of ', self.fullName)
        f_skeleton = open('{0}/defects/{1}/skeletons/skeletons_{2}.txt'.format(self.savePath, self.name, self.fullName), 'r')
        f_skin = open('{0}/defects/{1}/bodies/skins_{2}.txt'.format(self.savePath, self.name, self.fullName), 'r')

        count = 0
        skeleton = []
        skeletons = []
        for line in f_skeleton.readlines():
            line = line.strip().split(' ')
            # print(line)
            count += len(line)
            skeleton.append([int(line[j]) for j in range(len(line))])
            for elt in line:
                skeletons.append(int(elt))
        skin = []
        skins = []
        for line in f_skin.readlines():
            line = line.strip().split(' ')
            count += len(line)
            skin.append([int(line[j]) for j in range(len(line))])
            for elt in line:
                skins.append(int(elt))
        assert len(skeleton) == len(skin)

        body = []
        for i in range(len(skin)):
            for defect in skeleton[i]:
                body.append(defect)
            for deformed_hexagon in skin[i]:
                if deformed_hexagon not in body:
                    body.append(deformed_hexagon)

        nlist_body = defaultdict(list)
        for n_i, nb in enumerate(self.voro.nlist):
            ii, jj = nb
            if ii in body and jj in body:
                nlist_body[ii].append(jj)

        already_grouped = []
        agglomerates = []
        # for b1 in skeletons:
        for j in range(len(skeleton)):
            cj = []
            for b1 in skeleton[j]:
                if b1 not in already_grouped:
                    cj.append(b1)
                    already_grouped.append(b1)
                    for b2 in nlist_body[b1]:
                        if b2 not in already_grouped:
                            cj.append(b2)
                        already_grouped.append(b2)
            already_grouped = list(set(already_grouped + cj))
            for b3 in cj:
                for n3 in nlist_body[b3]:
                    if n3 in body and n3 not in already_grouped:
                        cj.append(n3)
                        already_grouped.append(n3)
            already_grouped = list(set(already_grouped + cj))
            if cj:
                agglomerates.append(cj)

        agglomerates = reversed(sorted(agglomerates, key=len))
        f_bodies = open('{0}/defects/{1}/bodies/bodies_{2}.txt'.format(self.savePath, self.name, self.fullName), 'w')
        for connected_component in agglomerates:
            for elt in connected_component:
                f_bodies.write('{} '.format(elt))
            f_bodies.write('\n')
        f_bodies.close()


class StructureAnalysis(object):

    def __init__(self, coord, L_x, L_y, min_k, max_k, num_kbins, ppName, name, save_path):
        self.coord = coord
        self.numPoints = coord.shape[0]
        self.L_x = L_x
        self.L_y = L_y
        self.min_k = min_k
        self.max_k = max_k
        self.num_kbins = num_kbins
        self.fullName = name
        self.resultPath = save_path
        self.P_name = ppName

    def compute_structure_factor_for_2d_data(self, file_name, boxsz=None, save=True):
        """ This function is defined by Michael A. Klatt (https://mklatt.org/) """
        
        system('mkdir {}/SF'.format(self.resultPath))
        num_points = self.coord.shape[0]
        if boxsz == None:
            boxsz = np.sqrt(num_points)
        else:
            boxsz = boxsz
        max_N = int((self.max_k - self.min_k) / (2 * np.pi / boxsz)) + 1
        N_k = (2 * max_N) ** 2
        k = np.zeros((N_k))
        kx = np.zeros((N_k))
        ky = np.zeros((N_k))
        nxny = np.zeros((N_k, 2))
        delta_k_0 = 2 * np.pi / boxsz
        delta_k_1 = 2 * np.pi / boxsz

        count = 0
        tmp = np.zeros((2, 1))
        for nx in range(1, max_N):
            ny = 0
            tmp[0] = nx * delta_k_0
            tmp[1] = ny * delta_k_1
            tmp_k = np.linalg.norm(tmp)
            if tmp_k < self.max_k and tmp_k > self.min_k:
                kx[count] = tmp[0]
                ky[count] = tmp[1]
                k[count] = tmp_k
                nxny[count, 0] = nx
                nxny[count, 1] = ny
                count += 1

        for nx in range(-max_N, max_N):
            for ny in range(1, max_N):
                tmp[0] = nx * delta_k_0
                tmp[1] = ny * delta_k_1
                tmp_k = np.linalg.norm(tmp)
                if tmp_k < self.max_k and tmp_k > self.min_k:
                    kx[count] = tmp[0]
                    ky[count] = tmp[1]
                    k[count] = tmp_k
                    nxny[count, 0] = nx
                    nxny[count, 1] = ny
                    count += 1
        kx = kx[:count]
        ky = ky[:count]
        k = k[:count]
        nxny = nxny[:count, :]

        # print('Determining structure factor ...')
        S = np.zeros((count), dtype=complex)  # structure factor
        # if self.tqdm_available:
        # for sj in tqdm(range(num_points)):
        #     S += np.exp(-1j * (kx * float(self.coord[sj, 0]) + ky * float(self.coord[sj, 1])))
        # else:
        for sj in range(num_points):
            S += np.exp(-1j * (kx * float(self.coord[sj, 0]) + ky * float(self.coord[sj, 1])))
            if sj % 1000 == 0:
                print('Computing SF for {} : Point {}'.format(self.fullName, sj))

        # The final (unbinned) structure factor for each absolute value of k
        unbinned = np.zeros((count, 4))
        unbinned[:, 0] = k
        unbinned[:, 1] = np.square(np.absolute(S)) / num_points
        unbinned[:, 2:] = nxny

        # # Output scattering intensity?
        # if args.sctint:
        #     np.savetxt(coords_file[:-4] + ".sctint", unbinned, fmt='%.10f %.10f %i %i %i')

        # Binning to estimate structure factor
        # if log_binning:
        #     if min_k == 0.0:
        #         print("Logarithmic binning automatically uses minimum wavenumber 0.01 instead of 0.")
        #         min_k = 0.01
        #     b = 10.0
        #     bins = np.logspace(np.log(min_k) / np.log(b), np.log(max_k) / np.log(b), num=num_kbins + 1, base=b)
        # else:
        bins = np.linspace(self.min_k, self.max_k + 1, self.num_kbins + 1)

        # binwidth = bins[1:] - bins[:-1]
        idx = np.digitize(unbinned[:, 0], bins)

        # output should contain
        # k(mean of bin)   <S>   err(S)   N_entries_in_this_bin
        output = np.zeros((np.size(bins) - 1, 4))
        output[:, 0] = (bins[1:] + bins[:-1]) * 0.5
        for i, bin_j in enumerate(idx):
            # i is row in unbinned, bin_j is bin to which it belongs
            output[bin_j - 1, 1] += unbinned[i, 1]
            output[bin_j - 1, 2] += (unbinned[i, 1]) ** 2
            output[bin_j - 1, 3] += 1

        # only evaluate bins with more than one k-value
        idx = np.where(output[:, 3] > 1)[0]

        output[idx, 2] = np.sqrt(np.fabs(
            output[idx, 2] / (output[idx, 3] - 1) / output[idx, 3] - np.power(output[idx, 1], 2) / output[
                idx, 3] / (
                    output[idx, 3] - 1) / output[idx, 3]))
        output[idx, 1] /= output[idx, 3]

        # for all others assign nan
        idx = np.where(output[:, 3] < 2)[0]
        output[idx, 1:3] = np.nan

        # Output structure factor
        if save:
            np.savetxt("{0}/SF/{1}/{2}.dat".format(self.resultPath, self.P_name, file_name), output, fmt="%.10f %.10f %.10f %i")
            # np.savetxt(coords_file[:-4] + "-Sk.dat", output, fmt="%.10f %.10f %.10f %i")
        else:
            return np.array(output)


    def structure_factor_for_3d_data(self, num_points, min_k, max_k, num_kbins, box_size=None):
        """ This function is defined by Michael A. Klatt (https://mklatt.org/) """
        
        if box_size == None:
            box_size = num_points ** (1. / 3)

        max_N = int((max_k - min_k) / (2 * np.pi / box_size)) + 1

        N_k = (2 * max_N) ** 3

        k = np.zeros((N_k))
        kx = np.zeros((N_k))
        ky = np.zeros((N_k))
        kz = np.zeros((N_k))
        nxnynz = np.zeros((N_k, 3))

        delta_k_0 = 2 * np.pi / box_size
        delta_k_1 = 2 * np.pi / box_size
        delta_k_2 = 2 * np.pi / box_size

        count = 0
        tmp = np.zeros((3, 1))
        for nx in range(-max_N, max_N):
            for ny in range(-max_N, max_N):
                nz = 0

                if ny > -nx or (ny == -nx and nx > 0):
                    tmp[0] = nx * delta_k_0
                    tmp[1] = ny * delta_k_1
                    tmp[2] = nz * delta_k_2

                    tmp_k = np.linalg.norm(tmp)

                    if tmp_k < max_k and tmp_k > min_k:
                        kx[count] = tmp[0]
                        ky[count] = tmp[1]
                        kz[count] = tmp[2]
                        k[count] = tmp_k

                        nxnynz[count, 0] = nx
                        nxnynz[count, 1] = ny
                        nxnynz[count, 2] = nz

                        count += 1
        for nx in range(-max_N, max_N):
            for ny in range(-max_N, max_N):
                for nz in range(1, max_N):

                    tmp[0] = nx * delta_k_0
                    tmp[1] = ny * delta_k_1
                    tmp[2] = nz * delta_k_2

                    tmp_k = np.linalg.norm(tmp)

                    if tmp_k < max_k and tmp_k > min_k:
                        kx[count] = tmp[0]
                        ky[count] = tmp[1]
                        kz[count] = tmp[2]
                        k[count] = tmp_k

                        nxnynz[count, 0] = nx
                        nxnynz[count, 1] = ny
                        nxnynz[count, 2] = nz

                        count += 1

        kx = kx[:count]
        ky = ky[:count]
        kz = kz[:count]
        k = k[:count]
        nxnynz = nxnynz[:count, :]

        # print('Determining structure factor ...')
        S = np.zeros((count), dtype=complex)  # structure factor
        if self.tqdm_available:
            for sj in tqdm(range(num_points)):
                S += np.exp(
                    -1j * (kx * float(self.coord[sj, 0]) + ky * float(self.coord[sj, 1]) + kz * float(
                        self.coord[sj, 2])))
        else:
            for sj in range(num_points):
                S += np.exp(-1j * (kx * self.coord[sj, 0] + ky * self.coord[sj, 1] + kz * self.coord[sj, 2]))

        # The final (unbinned) structure factor for each absolute value of k
        unbinned = np.zeros((count, 5))
        unbinned[:, 0] = k
        unbinned[:, 1] = np.square(np.absolute(S)) / num_points
        unbinned[:, 2:] = nxnynz

        # # Output scattering intensity?
        # if args.sctint:
        #     np.savetxt(coords_file[:-4] + ".sctint", unbinned, fmt='%.10f %.10f %i %i %i')

        # Binning to estimate structure factor
        # if log_binning:
        #     if min_k == 0.0:
        #         print("Logarithmic binning automatically uses minimum wavenumber 0.01 instead of 0.")
        #         min_k = 0.01
        #     b = 10.0
        #     bins = np.logspace(np.log(min_k) / np.log(b), np.log(max_k) / np.log(b), num=num_kbins + 1, base=b)
        # else:
        bins = np.linspace(min_k, max_k + 1, num_kbins + 1)

        # binwidth = bins[1:] - bins[:-1]
        idx = np.digitize(unbinned[:, 0], bins)

        # output should contain
        # k(mean of bin)   <S>   err(S)   N_entries_in_this_bin
        output = np.zeros((np.size(bins) - 1, 4))
        output[:, 0] = (bins[1:] + bins[:-1]) * 0.5
        for i, bin_j in enumerate(idx):
            # i is row in unbinned, bin_j is bin to which it belongs
            output[bin_j - 1, 1] += unbinned[i, 1]
            output[bin_j - 1, 2] += (unbinned[i, 1]) ** 2
            output[bin_j - 1, 3] += 1

        # only evaluate bins with more than one k-value
        idx = np.where(output[:, 3] > 1)[0]

        output[idx, 2] = np.sqrt(np.fabs(
            output[idx, 2] / (output[idx, 3] - 1) / output[idx, 3] - np.power(output[idx, 1], 2) / output[
                idx, 3] / (
                    output[idx, 3] - 1) / output[idx, 3]))
        output[idx, 1] /= output[idx, 3]

        # for all others assign nan
        idx = np.where(output[:, 3] < 2)[0]
        output[idx, 1:3] = np.nan

        # Output structure factor
        np.savetxt("%s/SF/%s/SF_q%.1f_%ibins_" % (self.resultPath, self.P_name, max_k, num_kbins) + name + ".dat", output,
                   fmt="%.10f %.10f %.10f %i")
        # np.savetxt(coords_file[:-4] + "-Sk.dat", output, fmt="%.10f %.10f %.10f %i")
        # return output


