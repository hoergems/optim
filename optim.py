import numpy as np
from sympy import *
from sympy.physics.quantum import TensorProduct

class Optim:
    def __init__(self):
        waypoints = self.get_waypoints()
        bad_points = self.get_bad_points()
        K = self.K(waypoints)        
        f, f2, epsi, symb_arr, A = self.get_f(waypoints, K)        
        self.optimize(waypoints, f[0], A, epsi, symb_arr, f2)        
        
    def optimize(self, waypoints, f, A, epsi, symb_arr, f2=None):       
        vals = Matrix(flatten([waypoints[i] for i in xrange(1, len(waypoints) - 1)]))        
        A_inv = A.inv() 
        gradient = self.gradient(f, symb_arr)
        delta = 10000.0
        ev = 10000.0
        while true:
            grad = gradient            
            for j in xrange(len(symb_arr)):
                print symb_arr
                print vals
                print flatten(epsi)
                grad = grad.subs(symb_arr[j], vals[j]) 
                       
            vals = vals - (1.0 / 2.0) * A_inv * grad
            
            ev = self.eval_f(f, epsi_symb_vec, vals)
            print ev
            if not f2 == None:
                f2_temp = f2
                for i in xrange(len(symb_arr)):
                    f2_temp = f2_temp.subs(symb_arr[i], vals[i])
                print "ev2 " + str(f2_temp)
                print "============================"
            if ev < 10e-8:
                break        
        idx = 0
        wp = []
        for i in xrange(0, len(vals), len(waypoints[0])):
            wp.append(np.array([float(vals[i + j]) for j in xrange(len(waypoints[0]))]))
        print wp
        sleep
        print vals
    
    def get_waypoints(self):
        x0 = np.array([0.0, 0.0])
        x1 = np.array([10.0, 10.0])        
        waypoints = [x0]        
        for i in xrange(2):
            s = np.random.uniform(0, 10, 2)            
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
        K = Matrix([[1, 0, 0, 0],
                    [-1, 1, 0, 0],
                    [0, 0, -1, 1],
                    [0, 0, 0, -1]])
        
        return K
        #I = Matrix.eye(len(waypoints[0]))
        #K_t = TensorProduct(K, I)
        print K_t.shape
        sleep
        K = Matrix.eye(len(waypoints[0]))
        
        print K
        sleep      
        K = Matrix.zeros(len(waypoints) - 2, len(waypoints) - 2)
        K[0, 0] = 1.0
        K[K.rows - 1, K.cols - 1] = -1.0
        for i in xrange(1, K.rows - 1):
            K[i, i - 1] = -1.0
            K[i, i] = 1.0
        return K
    
    def eval_f(self, f, epsi_symb_vec, epsi_val_vec):                
        f = lambdify(epsi_symb_vec, f, 'numpy')        
        return f(*flatten(epsi_val_vec))
        
    def get_f(self, waypoints, K):
        I = Matrix.eye(len(waypoints[0]))
        K = TensorProduct(K, I)
        epsi = [] 
        epsi_symb_vec = [Symbol("x_0_" + str(i)) for i in xrange(len(waypoints[0]))]
        for i in xrange(1, len(waypoints) - 1):
            for j in xrange(len(waypoints[0])):                
                epsi.append(Symbol("x_" + str(i) + "_" + str(j)))
                epsi_symb_vec.append(epsi[-1])                                 
        epsi_vec = Matrix(2 * epsi)
        epsi_symb_vec.extend([Symbol("x_" + str(len(waypoints) - 1) + "_" + str(i)) for i in xrange(len(waypoints[0]))])
        e = [-Symbol("x_0_" + str(i)) for i in xrange(len(waypoints[0]))]            
        #e = [-waypoints[0][i] for i in xrange(len(waypoints[0]))]
        for i in xrange(len(waypoints) - 2):
            for j in xrange(len(waypoints[i])):
                e.append(0.0)
        e.extend([Symbol("x_" + str(len(waypoints) - 1) + "_" + str(j)) for j in xrange(len(waypoints[0]))])        
        #e.extend([waypoints[len(waypoints) - 1][i] for i in xrange(len(waypoints[0]))])                 
        e = Matrix(flatten(e))
        
        f1 = 0.5 * epsi_vec.T * K.T * K * epsi_vec + epsi_vec.T * K.T * e + 0.5 * e.T * e
        sum = 0.0        
        epsi_vec2 = [[Symbol("x_" + str(i) + "_" + str(j)) for j in xrange(len(waypoints[0]))] for i in xrange(len(waypoints))]
        
        
        for i in xrange(1, len(epsi_vec2)):
            s1 = 0.0
            for j in xrange(len(waypoints[0])):
                s1 += (epsi_vec2[i][j] - epsi_vec2[i - 1][j])**2
            sum += s1
        f2 = (1.0 / 2.0) * sum                       
        return f1, f2, epsi_vec, epsi_symb_vec, K.T * K
    
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