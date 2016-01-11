import numpy as np
from sympy import *
from sympy.physics.quantum import TensorProduct

class Optim:
    def __init__(self):
        waypoints = self.get_waypoints()
        bad_points = self.get_bad_points()
        K = self.K(waypoints)
        
        f, f2, epsi_symb_vec = self.get_f(waypoints, K)
        self.eval_f(f, f2, epsi_symb_vec, waypoints)
        print self.K(waypoints)
    
    def get_waypoints(self):
        x0 = np.array([0.0, 0.0])
        x1 = np.array([10.0, 10.0])
        waypoints = [x0]
        x = np.array([0.0, 0.0])
        for i in xrange(3):
            x = x + np.array([1.0, 1.0])
            waypoints.append(x)
        return waypoints
    
    def get_bad_points(self):
        bp = [np.array([2.0, 3.0]),
              np.array([8.0, 7.0])]
        return bp
    
    def K(self, waypoints):        
        K = Matrix.zeros(len(waypoints), len(waypoints))        
        K[0, 0] = -1.0
        K[K.rows - 1, K.cols - 1] = -1.0
        for i in xrange(1, K.rows - 1):
            K[i, i - 1] = -1.0
            K[i, i] = 1.0
        return K
    
    def eval_f(self, f, f2, epsi_symb_vec, epsi_val_vec):        
        f = lambdify(epsi_symb_vec, f[0, 0], 'numpy')        
        print f(*flatten(epsi_val_vec))
        print f2
        sleep
        k = f.subs(epsi_symb_vec[0], Matrix([epsi_val_vec[0][0],
                                             epsi_val_vec[0][1]]))
        print k
        for i in xrange(len(epsi_symb_vec)):
            pass
        sleep
        
    def get_f(self, waypoints, K):
        I = Matrix.eye(len(waypoints[0]))
        K_t = TensorProduct(K, I)               
        epsi = []
        for i in xrange(len(waypoints)):
            ep = []
            for j in xrange(len(waypoints[i])):                
                x_str = "x_" + str(i) + "_" + str(j)
                ep.append(Symbol(x_str))
            epsi.extend(ep)               
        epsi_vec = Matrix(epsi)     
        
        A = K_t.T * K_t
        e = [[0.0 for i in xrange(len(waypoints[0]))] for j in xrange(len(waypoints))]
        e[0] = [-waypoints[0][i] for i in xrange(len(waypoints[0]))]
        e[len(e) - 1] = [waypoints[len(e) - 1][i] for i in xrange(len(waypoints[len(e) - 1]))]
        print waypoints
        print e
        e = Matrix(flatten(e))        
        b = K_t.T * e
        c = (1.0 / 2.0) * e.T * e 
        print c   
        f1 = (1.0 / 2.0) * epsi_vec.T * A * epsi_vec + epsi_vec.T * b + c
        
        sum = 0.0
        for i in xrange(1, len(waypoints)):
            s1 = 0.0
            for j in xrange(len(waypoints[i])):
                s1 += np.square(waypoints[i][j] - waypoints[i - 1][j])
            sum += s1
        f2 = (1.0 / 2.0) * sum        
        return f1, f2, epsi
    
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