import numpy as np
from sympy import *
from sympy.physics.quantum import TensorProduct

class Optim:
    def __init__(self):
        waypoints = self.get_waypoints()
        bad_points = self.get_bad_points()
        K = self.K(waypoints)        
        f, f2, epsi_symb_vec, A = self.get_f(waypoints, K)               
        symb_arr = flatten([[epsi_symb_vec[i][j] for j in xrange(len(epsi_symb_vec[i]))] for i in xrange(len(epsi_symb_vec))])              
        self.optimize(waypoints, f[0], symb_arr, A, epsi_symb_vec, f2)        
        
    def optimize(self, waypoints, f, symb_arr, A, epsi_symb_vec, f2=None):
        print f
        sleep
        vals = Matrix(flatten([waypoints[i] for i in xrange(1, len(waypoints) - 1)]))        
        A_inv = A.inv()  
        #print f.expand()
        #f2 = f2.subs(symb_arr[len(symb_arr) - 2], 10.0)
        #f2 = f2.subs(symb_arr[len(symb_arr) - 1], 10.0)
        #print " "
        #print f2.expand()
        #sleep  
        gradient = self.gradient(f, symb_arr)
        delta = 10000.0
        ev = 10000.0
        while true:
            grad = gradient            
            for j in xrange(len(symb_arr)):
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
        K = Matrix.zeros(len(waypoints) - 2, len(waypoints) - 2)
        K[0, 0] = 1.0
        K[K.rows - 1, K.cols - 1] = -1.0
        for i in xrange(1, K.rows - 1):
            K[i, i - 1] = -1.0
            K[i, i] = 1.0 
        print K                 
        return K
    
    def eval_f(self, f, epsi_symb_vec, epsi_val_vec):                
        f = lambdify(epsi_symb_vec, f, 'numpy')        
        return f(*flatten(epsi_val_vec))
        
    def get_f(self, waypoints, K):
        I = Matrix.eye(len(waypoints[0]))                
        K_t = TensorProduct(K, I)            
        #K_t = K         
        epsi = []
        for i in xrange(1, len(waypoints) - 1):
            ep = []
            for j in xrange(len(waypoints[i])):                
                x_str = "x_" + str(i) + "_" + str(j)
                ep.append(Symbol(x_str))
            epsi.append(Matrix(ep))                          
        epsi_vec = Matrix(epsi)        
        A = K_t.T * K_t        
        e = [-waypoints[0][i] for i in xrange(len(waypoints[0]))]
        for i in xrange(len(waypoints) - 4):
            for j in xrange(len(waypoints[i])):
                e.append(0.0)
        e.extend([waypoints[len(waypoints) - 1][i] for i in xrange(len(waypoints[0]))])                 
        e = Matrix(flatten(e))               
        a = (1.0 / 2.0) * epsi_vec.T * A * epsi_vec
        b = K_t.T * e
        c = (1.0 / 2.0) * e.T * e      
        f1 = a + epsi_vec.T * b + c 
                     
        sum = 0.0        
        epsi_vec2 = [[Symbol("x_" + str(i) + "_" + str(j)) for j in xrange(len(waypoints[0]))] for i in xrange(len(waypoints))]
        
        for i in xrange(1, len(epsi_vec2)):
            s1 = 0.0
            for j in xrange(len(waypoints[0])):
                s1 += np.square(epsi_vec2[i][j] - epsi_vec2[i - 1][j])
            sum += s1
        f2 = (1.0/ 2.0) * sum
        
        for i in xrange(len(waypoints[0])):
            f2 = f2.subs(epsi_vec2[0][i], waypoints[0][i])
            f2 = f2.subs(epsi_vec2[len(epsi_vec2) - 1][i], waypoints[len(waypoints) - 1][i])                              
        return f1, f2, epsi, A
    
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