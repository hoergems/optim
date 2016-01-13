import numpy as np
from sympy import *
from sympy.physics.quantum import TensorProduct
from plot import plot_2d_n_sets

class Optim:
    def __init__(self):
        waypoints = self.get_waypoints(25)
        obstacles = self.get_bad_points()
        print "Get K"
        K = self.K(waypoints) 
        print "Got K        "       
        f_prior, epsi, A = self.get_f_prior(waypoints, K)
        f_obst, obst_symb_vec = self.get_f_obst(waypoints)
        
        print "Got f"        
        trajectory = self.optimize(waypoints, obstacles, f_prior, f_obst, A, epsi, obst_symb_vec)
        print trajectory
        print obstacles
        plot_sets = [np.array(trajectory), np.array(obstacles)]
        
        
        plot_2d_n_sets(np.array(plot_sets),
                       x_range=[0.0, 15.0],
                       y_range=[0.0, 15.0],
                       plot_type="points")
        
    def get_f_obst(self, waypoints):
        dim = len(waypoints[0])
        epsi = []
        obst = []
        for i in xrange(1, len(waypoints) - 1):
            for j in xrange(dim):
                epsi.append(Symbol("x_" + str(i) + "_" + str(j)))
                obst.append(Symbol("o_" + str(i) + "_" + str(j)))
        epsi_vec = Matrix(epsi)
        obst_vec = Matrix(obst)
        
        f = 0.5 * epsi_vec.T * epsi_vec - epsi_vec.T * obst_vec + 0.5 * obst_vec.T * obst_vec
        #f = epsi_vec.T * epsi_vec - 2 * epsi_vec.T * obst_vec + obst_vec.T * obst_vec
        f = 1 / f[0]        
        return f, obst_vec
    
    def get_closest_obstacle(self, x, obstacles):        
        closest_obstacle = obstacles[0]
        min_dist = np.linalg.norm(x - closest_obstacle)
        for i in xrange(1, len(obstacles)):
            dist = np.linalg.norm(x - obstacles[i])
            if dist < min_dist:
                closest_obstacle = obstacles[i]
                min_dist = dist
        return closest_obstacle
    
    def substitude_closest_obstacle(self, f, vals, obstacles, dim):        
        v = flatten(vals)
        xs = []
        for i in xrange(0, len(v), dim):
            x = []
            for j in xrange(dim):
                x.append(float(v[i + j]))
            x = np.array(x)
            xs.append(x)
        o_vals = []
        for i in xrange(len(xs)):
            o = self.get_closest_obstacle(xs[i], obstacles)
            o_vals.append(o)
            for j in xrange(len(o)):
                f = f.subs("o_" + str(i + 1) + "_" + str(j), o[j])
        return f, flatten(o_vals)   
        
        
        sleep
        closest_obstacle = obstacles[0]
        for i in xrange(1, len(obstacles)):
            pass                 
        
    def optimize(self, waypoints, obstacles, f_prior, f_obst, A, epsi, obst_symb_vec, f2=None):
        obst_symb_vec = flatten(obst_symb_vec) 
        dim = len(waypoints[0])            
        vals = Matrix(flatten([waypoints[i] for i in xrange(1, len(waypoints) - 1)]))
        symb_arr = flatten(epsi) 
        print "invert A"       
        A_inv = A.inv()
        print "Get gradient"
        f1 = 0.5 * f_prior
        f2 =1000.0 * f_obst 
        f = f1 + f2 
        #f = 20.0 * f_obst       
        gradient = self.gradient(f, symb_arr)        
        dist = 100000000000000.0
        ev = 1000000000000.0
        
        theta = 1.0
        while true:
            grad = gradient
            grad, o_vals = self.substitude_closest_obstacle(grad, vals, obstacles, dim)            
            for j in xrange(len(symb_arr)):
                grad = grad.subs(symb_arr[j], vals[j]) 
                       
            vals = vals - theta * A_inv * grad
            ev_old = ev
            ev_prior = self.eval_f(f1, symb_arr, vals, obst_symb_vec, o_vals)
            ev_obst = self.eval_f(f2, symb_arr, vals, obst_symb_vec, o_vals)
            ev = self.eval_f(f, symb_arr, vals, obst_symb_vec, o_vals)
            print "ev_prior " + str(ev_prior)
            print "ev_obst " + str(ev_obst)            
            print "ev " + str(ev)
            old_dist = dist
            dist = ev_old - ev
            print "dist " + str(dist)
            print "======================"
            if dist < 0 or dist > old_dist:
                theta *= 0.5
            elif dist < 10e-5:
                break        
        trajectory = [waypoints[0]]
        for i in xrange(0, len(vals), dim):
            tr = []
            for j in xrange(dim):
                tr.append(float(vals[i + j]))
            trajectory.append(np.array(tr))
        trajectory.append(waypoints[-1]) 
        return trajectory       
        
    
    def get_waypoints(self, num=2):
        x0 = np.array([0.0, 0.0])
        x1 = np.array([10.0, 10.0])        
        waypoints = [x0]        
        for i in xrange(num):
            s = np.random.uniform(0, 10, len(x0))            
            
            waypoints.append(s)
        waypoints.append(x1)        
        return waypoints
    
    def gradient(self, f, variables):
        """ Compute the gradient of a function        
        """
        
        grad = Matrix([diff(f, variables[i]) for i in xrange(len(variables))]) 
        return grad
    
    def get_bad_points(self):
        bp = [np.array([2.0, 2.0]),
              np.array([8.0, 7.0]),
              np.array([1.0, 8.0]),
              np.array([8.0, 1.0]),
              np.array([8.0, 14.0])]
        return bp
    
    def K(self, waypoints):
        d = len(waypoints[0])
        n = len(waypoints)
        K = Matrix.zeros(d * (n - 1), d * (n - 2))
        for i in xrange(K.rows - d):
            K[i, i] = 1        
        l = K.cols - 1
        for i in xrange(K.rows - 1, d - 1, -1):
            K[i, l] = -1
            l -= 1        
        return K
    
    def eval_f(self, f, epsi_symb_vec, epsi_val_vec, obst_symb_vec, obst_val_vec):           
        for i in xrange(len(epsi_symb_vec)):
            f = f.subs(epsi_symb_vec[i], epsi_val_vec[i])
        for i in xrange(len(obst_symb_vec)):
            f = f.subs(obst_symb_vec[i], obst_val_vec[i])
        return f
        
    def get_f_prior(self, waypoints, K):
        I = Matrix.eye(len(waypoints[0]))        
        epsi = []        
        for i in xrange(1, len(waypoints) - 1):
            for j in xrange(len(waypoints[0])):                
                epsi.append(Symbol("x_" + str(i) + "_" + str(j)))                                                        
        epsi_vec = Matrix(epsi)        
        e = [0 for i in xrange(K.shape[0])]
        for i in xrange(len(waypoints[0])):
            e[i] = -Symbol("x_0_" + str(i))
            e[len(e) - len(waypoints[0]) + i] = Symbol("x_" + str(len(waypoints) - 1) + "_" + str(i))             
        e = Matrix(flatten(e))
        f1 = 0.5 * epsi_vec.T * (K.T * K) * epsi_vec + epsi_vec.T * K.T * e + 0.5 * e.T * e        
        e_v = flatten(e)
        for i in xrange(len(waypoints[0])):            
            f1 = f1.subs(-e_v[i], waypoints[0][i])
            f1 = f1.subs(e_v[len(e_v) - len(waypoints[0]) + i], waypoints[len(waypoints) - 1][i])                          
        return f1[0], epsi_vec, K.T * K
    
    def calc_f_dist(self, waypoints, bad_points):
        for waypoint in waypoints:
            index, min_dist = self.get_nearest_bad_points(waypoint, bad_points)
            print index, min_dist
        
    def get_nearest_bad_points(self, waypoint, bad_points):
        index = 0
        min_dist = 1000000.0
        for i in xrange(len(bad_points)):
            dist = np.linalg.norm(bad_points[i] - waypoint)
            if dist < min_dist:
                min_dist = dist
                index = i
        return index, min_dist
        

if __name__ == "__main__":
    Optim()