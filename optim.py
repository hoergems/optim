import numpy as np
from sympy import *
from sympy.physics.quantum import TensorProduct

class Optim:
    def __init__(self):
        waypoints = self.get_waypoints(50)
        bad_points = self.get_bad_points()
        print "Get K"
        K = self.K(waypoints) 
        print "Got K        "       
        f, epsi, A = self.get_f(waypoints, K)
        print "Got f"
        self.optimize(waypoints, f[0], A, epsi)        
        
    def optimize(self, waypoints, f, A, epsi, f2=None):             
        vals = Matrix(flatten([waypoints[i] for i in xrange(1, len(waypoints) - 1)]))
        symb_arr = flatten(epsi) 
        print "invert A"       
        A_inv = A.inv()
        print "Get gradient" 
        gradient = self.gradient(f, symb_arr)
        print "got gradient"
        delta = 100000000000000.0
        ev = 1000000000000.0
        while true:
            grad = gradient            
            for j in xrange(len(symb_arr)):
                grad = grad.subs(symb_arr[j], vals[j])            
            vals = vals - 0.2 * A_inv * grad
            ev_old = ev
            ev = self.eval_f(f, symb_arr, vals)
            print "ev old " + str(ev_old)
            print "ev " + str(ev)
            dist = ev_old - ev
            print dist
            if dist < 10e-8:
                break        
        print vals
    
    def get_waypoints(self, num=2):
        x0 = np.array([0.0, 0.0])
        x1 = np.array([10.0, 10.0])        
        waypoints = [x0]        
        for i in xrange(num):
            s = np.random.uniform(0, 10, len(x0))            
            x0 = x0 + s
            waypoints.append(x0)
        waypoints.append(x1)        
        return waypoints
    
    def gradient(self, f, variables):
        """ Compute the gradient of a function        
        """
        
        grad = Matrix([diff(f, variables[i]) for i in xrange(len(variables))]) 
        return grad
    
    def get_bad_points(self):
        bp = [np.array([2.0, 3.0]),
              np.array([8.0, 7.0])]
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
    
    def eval_f(self, f, epsi_symb_vec, epsi_val_vec):        
        for i in xrange(len(epsi_symb_vec)):
            f = f.subs(epsi_symb_vec[i], epsi_val_vec[i])
        return f
        
    def get_f(self, waypoints, K):
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
        return f1, epsi_vec, K.T * K
    
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