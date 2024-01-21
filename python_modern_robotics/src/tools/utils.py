import sympy as sym
import numpy as np
import matplotlib.pyplot as plt

from numpy.linalg import inv
from typing import List
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform



s,c,alpha, l = sym.symbols("s,c,alpha,l")

def exp_to_rotmatrix(vec)->np.matrix:
    mat = sym.Matrix([[0,-vec[2], vec[1]],
                      [vec[2], 0, -vec[0]],
                      [-vec[1], vec[0], 0]])
    I = sym.eye(3)
    
    exp_mat = I + mat*s + (mat**2)*(1-c)
    
    vec_mat = sym.Matrix([[vec[0], vec[1], vec[2]]]).transpose()
    ww_mat = vec_mat*vec_mat.transpose()

    
    k = (I - exp_mat)*mat + alpha*ww_mat
    return k

def exp_to_Hmatrix(twist)->np.matrix:
    w = twist[3:]
    v = twist[:3]
    I = sym.eye(3)
    
    rot = exp_to_rotmatrix(w)
    cross_product = (sym.Matrix(w).cross(sym.Matrix(v)))
    
    trans = (I - rot)*cross_product

    
class robot_kinematics:
    """
    class gets the coordinate frames and joint_types and returns the symbolic 
    forward kinematics matrix
    
    ...
    Attributes
    ----------
    frames: np.matrix(4x4)
        input are the frames and the joint types 
    
    joint_types: np.matrix
        type of joint in the frame

    """
    
    def __init__(self, frames) -> None:
        self.frames = [i[0] for i in frames]
        self.joint_types = [i[1] for i in frames]
        
        #initializing the class variables 
        self.no_of_joints = len(self.frames)
        self.H0 = self.frames[-1]
        self.H_array = []
        self.exp_twist_array = []
        self.Twist_array = [0]*len(frames)
        self.Jacobian_array = [0]*len(frames)
        
        s,c,alpha, l = sym.symbols("s,c,alpha,l")
        self.sin = sym.symbols("s1:%d"%(self.no_of_joints + 1))
        self.cos = sym.symbols("c1:%d"%(self.no_of_joints + 1))
    
    def tildeOperator(self, vec:np.array)->np.matrix:
        """
        gets the vector and converts to matrix
        
        parameters
        ----------
        vec: np.array[x,x,x]
            the vector that to be converted to matrix
        """
        
        mat = sym.Matrix([[0,-vec[2], vec[1]],
                      [vec[2], 0, -vec[0]],
                      [-vec[1], vec[0], 0]])
        
        return mat
        
    def getExpToRotationMatrix(self, vec:np.array, joint_no = 0)->np.matrix:
        """
        gets the vector around which axis of rotation is happening and converts to rotational matrix
        
        parameters
        ----------
        vec: np.array[x,x,x]
            the vector that to be converted to matrix
        """
        mat = self.tildeOperator(vec)
        I = sym.eye(3)
        
        rodrigues_formula = I + mat*self.sin[joint_no] + (mat**2)*(1-self.cos[joint_no]) #rodrigues formula

        return rodrigues_formula

    def getTwistToHmatrix(self, twist:np.array, joint_no)->np.matrix:
        """
        gets the twist[1x6] anf convert it into Homogenous matrix
        parameters
        ----------
        twist: np.array[1x6]
            the vector is taken innput to get converted to H[4x4] Homogenous matrix
        """
        rot_vec = sym.Matrix(twist[:3])
        trans_vec = sym.Matrix(twist[3:])
        I = sym.eye(3)

        rot_mat = self.tildeOperator(rot_vec) #vector is converted to matrix using tilde operator
        rot_component = self.getExpToRotationMatrix(rot_vec, joint_no)
        
        a = (sym.eye(3) - rot_component)
        b = rot_vec.cross(trans_vec)
        c = trans_vec.dot(rot_vec)*rot_vec.T
        w_norm = 1/rot_vec.dot(rot_vec)
        trans_component = w_norm*(a*b + c.T)
        
        
        twist_rodriguez = sym.eye(4)
        twist_rodriguez[0:3,0:3] = rot_component
        twist_rodriguez[0:3,-1] = trans_component
        
        
        
        return twist_rodriguez
    
    def brockettFormula(self) -> np.matrix:
        """
        computes the forward kinematics and store it in array

        Args:
            joint_index (int): nth index of the joint
        """
        #helper function
        def exp_mulitplication(index):
            """ gets the how many time we should mulitply the e^T

            Args:
                index (int): nth joint
            """
            H = np.eye(4)
            if index != 0 and index != 1:
                for i in range(index-1):
                    H*=self.exp_twist_array[i]
            else:
                pass
                    
            return H
        
        for joint in range(self.no_of_joints):
            #from Homogenous transformation we know that T_{s,a}(a wrt space(s)),
            # T_{a,b} = T_{s,a}*T_{s,b} -=> T_{a,b} = T_{a,s}^(-1)*T_{s,b} 
            
            if joint == 0:
                H_0 = np.eye(4)
                twist = self.twistComputation( self.frames[0], self.joint_types[joint])
                print("joint :",joint, "    twist :",twist) 
                e_T01 = self.getTwistToHmatrix(twist, joint)
                self.exp_twist_array.append(e_T01)
                H_q = e_T01@H_0
                self.Twist_array[0] = twist
                self.H_array.append(H_q)
                
            else:
               
                H_0 = self.frames[joint]
                T_ab_local = np.linalg.inv(self.frames[joint-1])@self.frames[joint]
                
                twist = self.twistComputation( T_ab_local, self.joint_types[joint]) 
                print("joint :",joint, "    twist :",twist)               
                T0_i = self.Adjoint(self.H_array[joint-1])@twist
                e_T0i = self.getTwistToHmatrix(T0_i, joint)
                self.exp_twist_array.append(e_T0i)
                self.Twist_array[joint] = T0_i 
                H_q = exp_mulitplication(joint+1)@H_0
                self.H_array.append(H_q)
                
          
        
        return self.H_array # returns the forward kinematics of the robot
        
    
    def geometricJacobian(self, H):
        pass
        
        
    def combineRotTranstoHMatrix(self, rot: np.matrix, translate:np.matrix) -> np.matrix:
        H_mat = np.eye(4) #initialize the matrix to be Identity(4)
        H_mat[0:3,0:3] = rot
        H_mat[0:3,-1] = translate
        return H_mat
    
    def twistComputation(self, next_frame:np.matrix, joint_type:str)->np.array:
        """gets the H_i and H_{i+1}, to compute the Twist T^i_{i-1}

        Args:
            origin (np.matrix): first H matrix 
            next_frame (np.matrix): next H matris
            joint_type (str): type of the jpint

        Returns:
            np.array: returns the constant twist [1x6]
        """
            
        def getRAndUnitOmega(next_frame):
            r = np.array(next_frame[0:3,-1])
            omega = np.array(next_frame[0:3,-2])
            
            length = 1 # if input frames are orthonormal norm is not needed
  
            unit_omega = omega/length
            
            return unit_omega.T,r.T
        
        
        omega,r = getRAndUnitOmega(np.matrix(next_frame))
        if joint_type == "R":
            #ROTATION omega,r
            # omega = next_frame[:3,-2]
           
            v = np.cross(r,omega) # v = rxw
            T_r = np.array([omega.item(0),omega.item(1),omega.item(2),v.item(0), v.item(1), v.item(2)])
            return T_r
            
        if joint_type == "P":
            #PRISMATIC
            T_p = np.array([0,0,0,omega.item(0),omega.item(1),omega.item(2)])
            return T_r
            
            
        if joint_type == "S":
            #SCREW
            pitch = 1
            v = np.cross(r,omega) # v = rxw
            T_s = np.array([omega.item(0),omega.item(1),omega.item(2),v.item(0), v.item(1), v.item(2) + pitch])# have to look into it (not correct)
            return T_s
    def H_zer0(self, joint_index:int)-> np.matrix:
        """ gives the H[0] matrix for the nth joint from the coordinate frame input"""
        return self.frames[joint_index]
        
    def HomogenousMatrixComputation(self, joint_index):
        pass
    
    def Adjoint(self, H:np.matrix)-> np.matrix:
        """gets the H(q) matrix and returns the 6x6 matrix from reference frame to inertial frame

        Args:
            H (np.matrix): homogenous matrix to nth joint
        
        Returns:
            adj_H[np.matrix]: 6x6 matrix
        """
    
        R = H[:3,:3]
        p = np.array(H[-1,:3])[0]
        p_tilde = self.tildeOperator(p)
        

        Adj = sym.eye(6)*0
        Adj[:3,:3] = R
        Adj[:3,3:] = sym.eye(3)*0
        Adj[3:,:3] = p_tilde*R
        Adj[3:,3:] = R

        return Adj
    
    def plotting_robot(self,joint_angles=None):
        
        #initializing the figure
        fig = plt.figure()
        self.ax = fig.add_subplot(111, projection='3d')
        self.ax.set_aspect('equal')
        self.ax.set_xlim(-3,7)
        self.ax.set_ylim(-3,7)
        self.ax.set_zlim(-3,7)
        
        self.ax.set_xlabel('x')
        self.ax.set_ylabel('y')
        self.ax.set_zlabel('z')
        
        #helper function
        # Python code to merge dictionary
        
        def Merge(dict1, dict2):
            for i in dict2.keys():
                dict1[i]=dict2[i]
            return dict1
        
        def ang_subs(index,joint_angles):
            """gets the index and sustitutes the angle values

            Args:
                index (int): index of the joint
            """
            sin_dict = {}
            cos_dict = {}
            for i in range(index):
                sin_dict[self.sin[i]] = np.sin(np.radians(joint_angles[i]))
                cos_dict[self.cos[i]] = np.cos(np.radians(joint_angles[i]))
            
            return Merge(sin_dict, cos_dict)
        
        def draw_line(p1,p2):
            x1,y1,z1 = p1[0],p1[1],p1[2]
            x2,y2,z2 = p2[0],p2[1],p2[2]
            self.ax.plot([x1, x2], [y1,y2],[z1,z2],linewidth=2, color="limegreen")
            self.ax.scatter(x1,y1,z1, c='red', marker='o', s=75)
            self.ax.scatter(x2, y2, z2, c='red', marker='o', s=75)
            
        def draw_axis(joint):
            H = self.frames[joint]
           
            origin = self.frames[joint]@np.array([[0,0,0,1]]).T
          
            x1,y1,z1 = origin.item(0),origin.item(1),origin.item(2)
            x_pt = np.array([[1,0,0,1]])
            y_pt = np.array([[0,1,0,1]])
            z_pt = np.array([[0,0,1,1]])
            x_arrow_tip = H@x_pt.T
            y_arrow_tip = H@y_pt.T
            z_arrow_tip = H@z_pt.T
                       
            
            
            # # Here we create the arrows:
            arrow_prop_dict = dict(mutation_scale=20, arrowstyle='->', shrinkA=0, shrinkB=0)

            x_ = Arrow3D(x1, y1, z1, x_arrow_tip.item(0), x_arrow_tip.item(1), x_arrow_tip.item(2), **arrow_prop_dict, color='r')
            self.ax.add_artist(x_)
            y_ = Arrow3D(x1, y1, z1, y_arrow_tip.item(0), y_arrow_tip.item(1), y_arrow_tip.item(2), **arrow_prop_dict, color='g')
            self.ax.add_artist(y_)
            z_ = Arrow3D(x1, y1, z1, z_arrow_tip.item(0), z_arrow_tip.item(1), z_arrow_tip.item(2), **arrow_prop_dict, color='b')
            self.ax.add_artist(z_)

            # Give them a name:
            
            # self.ax.text(x1 + 0.1, 0, 0, r'$x$')
            # self.ax.text(0, y1 + 0.1, 0, r'$y$')
            # self.ax.text(0, 0, z1 + 0.1, r'$z$')
                    
        
        if joint_angles != None :
            origin = [0,0,0]
            robot_joint_coordinates = []
            robot_joint_coordinates.append(origin)
            for joint in range(len(self.H_array)):
                data_dict = ang_subs(joint, joint_angles) # give the substitued values of sin and cosine of joint angles
                H = self.H_array[joint].subs(data_dict)
                draw_axis(joint)
                last_column = H[0:3,-1]        
                x,y,z = last_column[0], last_column[1],last_column[2]
                robot_joint_coordinates.append([x,y,z])
                
            for pt in range(len(robot_joint_coordinates) -1):
                draw_line(robot_joint_coordinates[pt],robot_joint_coordinates[pt+1])
                
                
    
        else:
            print("plotting the origin configuration, enter the angles correctly to plot the configuration")
        
        
  
            
                
class Arrow3D(FancyArrowPatch):
    
    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._xyz = (x, y, z)
        self._dxdydz = (dx, dy, dz)

    def draw(self, renderer):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (dx, dy, dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)
        
    def do_3d_projection(self, renderer=None):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))

        return np.min(zs)      
            
class Twist:
    def __init__(self, frame1:np.matrix, next_frame:np.matrix, joint_type:str) -> None:
        self.origin = frame1
        self.frame = next_frame
        self.joint_type = joint_type
        

        
    
if __name__ == "__main__":
    l1,l2,l3 =2,2,2 #link length
    first_joint = [np.eye(4),"R"]
    second_joint = [np.matrix([[1,0,0,0],
                               [0,0,-1,0],
                               [0,1,0,l1],
                               [0,0,0,1]]),"R"]
    third_joint = [np.matrix([[0,0,1,0],
                              [0,1,0,0],
                              [-1,0,0,l1+ l2],
                              [0,0,0,1]]), "R"]
    fourth_joint = [np.matrix([[1,0,0,0],
                              [0,1,0,0],
                              [0,0,1,l1+ l2 +l3],
                              [0,0,0,1]]), "R"]
    
    frames = [first_joint, second_joint, third_joint, fourth_joint]
    
    robot = robot_kinematics(frames)
    fk = robot.brockettFormula()
    print(fk)
    robot.plotting_robot([45,45,45,30])
    plt.show()
    
    
    
    


