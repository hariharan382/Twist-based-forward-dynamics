import matplotlib.pyplot as plt
import numpy as np

#library 
from tools.utils import robot_kinematics

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
    robot.plotting_robot([45,45,45,30])
    plt.show()