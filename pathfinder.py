import numpy as np
from math import *
from tkinter import *
from time import *
import threading
from scipy import special

EXPANSION_PRECISION = 15
MIN_RADIUS = 0.2
DENSITY_FACTOR = 1.8 # must be >1

ROBOT_SPEED = 0.5

WALLS_REPULSION = 1

SCREEN_SIZE = 800

UPDATE_FACTOR = 1.5

FRAMERATE = 10

class Obstacle():
    def __init__(self,pos_a,pos_b):
        self.pos_a = pos_a
        self.pos_b = pos_b
    def get_distance(self,pos):
        xp,yp = pos
        xa,ya = self.pos_a
        xb,yb = self.pos_b
        dist_ab = sqrt((xa-xb)**2+(ya-yb)**2)
        scal = ((xp-xa)*(xb-xa)+(yp-ya)*(yb-ya))/(dist_ab**2)
        if scal<0:
            return(sqrt((xa-xp)**2+(ya-yp)**2))
        elif scal>1:
            return(sqrt((xb-xp)**2+(yb-yp)**2))
        else:
            x = xa+scal*(xb-xa)
            y = ya+scal*(yb-ya)
            return(sqrt((x-xp)**2+(y-yp)**2))

class MeshPoint():
    def __init__(self,pos,world,custom_radius = 0):
        self.pos = pos
        self.world = world
        if custom_radius == 0:
            self.radius = self.get_closer_obstacle()/DENSITY_FACTOR
            self.valid = (DENSITY_FACTOR*self.radius >= MIN_RADIUS)
            self.radius = max(self.radius,MIN_RADIUS)
        else:
            self.radius = custom_radius
            self.valid = 0
        self.neighbors = []
        self.corrupted = 0
        self.djikstra_dist = float("inf")
        self.djikstra_previous = None
        self.djikstra_validated = 0
        self.neighbor_to_corrupted = 0
    def store_in_box(self):
        x,y = self.pos
        self.world.boxes[floor(x)][floor(y)].append(self)
    def get_closer_obstacle(self):
        min = float("inf")
        for obstacle in self.world.obstacles:
            dist = obstacle.get_distance(self.pos)
            if dist<min:
                min = dist
        return(min)
    def get_distance(self,point):
        x,y = self.pos
        xp,yp = point.pos
        return(sqrt((x-xp)**2+(y-yp)**2))
    def get_djikstra_distance(self,point):
        return(self.get_distance(point)+WALLS_REPULSION/(self.radius+point.radius))
    def check_valid(self,father):
        if self.valid:
            x,y = self.pos
            for point in father.neighbors:
                if self.get_distance(point)<point.radius:
                    return(0)
            return(1)
        else:
            return(0)
    def find_neighbors(self):
        result = []
        d = 1+floor(self.world.biggest_radius)
        # print(d)
        xs = floor(self.pos[0])
        ys = floor(self.pos[1])
        potential_neighbors = []
        for x in range(xs-d,xs+d+1):
            for y in range(ys-d,ys+d+1):
                if ((x>=0) and (x<10)) and ((y>=0) and (y<10)):
                    to_add = self.world.boxes[x][y]
                    for point in to_add:
                        potential_neighbors.append(point)
        for point in potential_neighbors:
            # print("1 point")
            if self.get_distance(point)<(self.radius+point.radius):
                # print(self.get_distance(point))
                # print("adding neighbor")
                self.neighbors.append(point)
    def link_to_neighbors(self):
        d = 1+floor(self.world.biggest_radius)
        # print(d)
        xs = floor(self.pos[0])
        ys = floor(self.pos[1])
        potential_neighbors = []
        for x in range(xs-d,xs+d+1):
            for y in range(ys-d,ys+d+1):
                if ((x>=0) and (x<10)) and ((y>=0) and (y<10)):
                    to_add = self.world.boxes[x][y]
                    for point in to_add:
                        potential_neighbors.append(point)
        for point in potential_neighbors:
            # print("1 point")
            if self.get_distance(point)<(self.radius+point.radius):
                # print(self.get_distance(point))
                # print("adding neighbor")
                self.neighbors.append(point)
                point.neighbors.append(self)
    def expand(self,points_to_expand):
        # print(self.pos)
        x,y = self.pos
        r = self.radius
        potential_new_points = []
        for teta in np.linspace(0,2*pi,EXPANSION_PRECISION):
            new_point = MeshPoint([x+r*1.001*cos(teta),y+r*1.001*sin(teta)],self.world,1)
            valid = 1
            for neighbor in self.neighbors:
                if neighbor.get_distance(new_point) < neighbor.radius:
                    valid = 0
            if valid:
                potential_new_points.append(MeshPoint([x+r*1.001*cos(teta),y+r*1.001*sin(teta)],self.world))
        potential_new_points_valid = []
        for point in potential_new_points:
            if point.check_valid(self):
                potential_new_points_valid.append(point)
            else: 
                ()
                # print("WRONG")
        # print(len(potential_new_points_valid),'efaeef')
        radiuses = []
        for point in potential_new_points_valid:
            radiuses.append(point.radius)
        order = np.argsort(-np.array(radiuses))
        sorted_new_points = []
        for o in order:
            sorted_new_points.append(potential_new_points_valid[o])
        count = 0
        # print("bbb")
        while(len(sorted_new_points)>0):
            # print("ccc")
            new_sorted_new_points = []
            chosen_point = sorted_new_points[0]
            
            if chosen_point.radius > self.world.biggest_radius:
                self.world.biggest_radius = chosen_point.radius
            chosen_point.link_to_neighbors()
            chosen_point.store_in_box()
            count += 1
            self.world.mesh.append(chosen_point)
            points_to_expand.append(chosen_point)
            if len(sorted_new_points)>1:
                for i in range(1,len(sorted_new_points)-1):
                    point = sorted_new_points[i]
                    if point.get_distance(chosen_point)>chosen_point.radius:
                        new_sorted_new_points.append(point)
            sorted_new_points = new_sorted_new_points
        return(count)

class Robot():
    def __init__(self,world):
        self.world = world
        world.robot = self
        self.pos = world.start.pos
        self.speed = ROBOT_SPEED
        self.goal = self.world.djikstra_path[-2]
        x,y = np.array(self.pos)*SCREEN_SIZE/10
        y = SCREEN_SIZE - y
        r = MIN_RADIUS*SCREEN_SIZE/10
        self.drawing = self.world.canvas.create_oval(x-r,y-r,x+r,y+r,fill = "red")
        # self.goal_drawing = self.world.canvas.create_oval(x-r,y-r,x+r,y+r,fill = "green")
        self.last_update = time()
        # self.update()
    def reset_goal(self):
        if len(self.world.djikstra_path)>1:
            self.goal = self.world.djikstra_path[-2]
        else:
            self.goal = None
    def update(self):
        dt = time()-self.last_update
        self.last_update = time()
        if self.goal:
            x,y = self.pos
            xg,yg = self.goal.pos
            if sqrt((x-xg)**2+(y-yg)**2)<self.speed*dt:
                self.world.start = self.goal
                self.world.djikstra_path.pop(-1)
                self.world.display_djikstra()
                self.reset_goal()
            else:
                dist = sqrt((xg-x)**2+(yg-y)**2)
                dx = ((xg-x)/dist)*self.speed*dt
                dy = ((yg-y)/dist)*self.speed*dt
                self.pos = [x+dx,y+dy]
        self.display()
        self.world.canvas.after(int(1000/FRAMERATE),self.update)
    def display(self):
        x,y = np.array(self.pos)*SCREEN_SIZE/10
        y = SCREEN_SIZE - y
        r = MIN_RADIUS*SCREEN_SIZE/10
        self.world.canvas.delete(self.drawing)
        self.drawing = self.world.canvas.create_oval(x-r,y-r,x+r,y+r,fill = "red")
        # if self.goal:
        #     x,y = np.array(self.goal.pos)*SCREEN_SIZE/10
        #     y = SCREEN_SIZE - y
        #     r = MIN_RADIUS*SCREEN_SIZE/10
        #     self.world.canvas.delete(self.goal_drawing)
        #     self.goal_drawing = self.world.canvas.create_oval(x-r,y-r,x+r,y+r,fill = "green")
        
    

class World():
    def __init__(self,obstacles,start,end):
        self.obstacles = obstacles
        self.start = MeshPoint(start,self)
        self.start.radius = MIN_RADIUS
        self.end = MeshPoint(end,self)
        self.end.radius = MIN_RADIUS
        self.mesh = [self.start,self.end]
        self.window = Tk()
        self.canvas = Canvas(self.window, width=SCREEN_SIZE, height=SCREEN_SIZE)
        self.canvas.focus_set()
        self.canvas.bind("<Button-1>", self.add_wall)
        self.canvas.bind("<Motion>", self.update_wall)
        self.canvas.bind("<ButtonRelease>", self.end_wall)
        self.canvas.pack()
        self.new_wall_start = []
        self.new_wall = None
        self.lines = []
        self.djikstra_lines = []
        self.djikstra_path = []
        self.boxes = [[[] for i in range(10)] for j in range(10) ]
        self.biggest_radius = MIN_RADIUS
        self.generate_mesh()
        self.original_mesh = [point.pos for point in self.mesh][2:]
        self.last_mesh_update = time()
        self.beziers = []
    def add_obstacle(self,event):
        self.new_wall_start = 1
        x,y = [event.x*10/SCREEN_SIZE,10-event.y*10/SCREEN_SIZE]
        self.new_wall = Obstacle([x,y],[x+0.1,y+0.1])
        self.obstacles.append(self.new_wall)
    def add_wall(self,event):
        x,y = [event.x*10/SCREEN_SIZE,10-event.y*10/SCREEN_SIZE]
        self.new_wall_start = [x,y]
    def update_obstacle(self,event):
        if self.new_wall_start:
            x,y = [event.x*10/SCREEN_SIZE,10-event.y*10/SCREEN_SIZE]
            self.new_wall.pos_a = [x,y]
            self.new_wall.pos_b = [x+0.1,y+0.1]
            self.start = MeshPoint(self.robot.pos,self)
            self.update_mesh()
            self.robot.reset_goal()
            self.display()
    def update_wall(self,event):
        if self.new_wall_start and time()-self.last_mesh_update>1/FRAMERATE:
            xs,ys = self.new_wall_start 
            x,y = [event.x*10/SCREEN_SIZE,10-event.y*10/SCREEN_SIZE]
            dist = sqrt((xs-x)**2+(ys-y)**2)
            if dist > 0.01:
                if self.new_wall:
                    self.new_wall.pos_b = [x,y]
                else:
                    self.new_wall = Obstacle(self.new_wall_start,[x,y])
                    self.obstacles.append(self.new_wall)
                self.start = MeshPoint(self.robot.pos,self)
                self.update_mesh()
                self.robot.reset_goal()
                self.display()
            self.last_mesh_update = time()
    def end_obstacle(self,event):
        self.new_wall = None
        self.new_wall_start = 0
        self.obstacles.pop(-1)
        self.start = MeshPoint(self.robot.pos,self)
        self.update_mesh()
        self.robot.reset_goal()
        self.display()
    def end_wall(self,event):
        self.new_wall = None
        self.new_wall_start = []
    def generate_mesh(self):
        # print("b")
        self.boxes = [[[] for i in range(10)] for j in range(10) ]
        self.biggest_radius = MIN_RADIUS
        t0 = time()
        points_to_expand = [self.start,self.end]
        self.start.store_in_box()
        self.end.store_in_box()
        n_points_to_expend = 2
        while (n_points_to_expend > 0):
            # print(n_points_to_expend)
            dn = points_to_expand[0].expand(points_to_expand)
            n_points_to_expend += dn 
            points_to_expand.pop(0)
            n_points_to_expend -= 1
        t1 = time()
        self.djikstra()
        print("mesh : ",t1-t0)
        print("djikstra : ",time()-t1)
        print("---")
    def update_mesh(self):
        # self.start.neighbors = []
        # self.end.neighbors = []
        # self.mesh = [self.start,self.end]
        # self.generate_mesh()
        
        self.boxes = [[[] for i in range(10)] for j in range(10) ]
        self.mesh = []
        self.start = MeshPoint(self.start.pos,self)
        self.end = MeshPoint(self.end.pos,self)
        self.start.store_in_box()
        self.end.store_in_box()
        self.start.radius = MIN_RADIUS
        self.end.radius = MIN_RADIUS
        self.mesh = [self.start,self.end]
        corrupted_points = []
        self.biggest_radius = MIN_RADIUS
        for coords in self.original_mesh:
            point = MeshPoint(coords,self)
            if point.get_closer_obstacle()<UPDATE_FACTOR*point.radius:
                corrupted_points.append(point)
                point.find_neighbors()
                point.corrupted = 1
            else:
                if point.radius > self.biggest_radius:
                    self.biggest_radius = point.radius
                self.mesh.append(point)
                point.store_in_box()
                point.link_to_neighbors()
        points_to_expand = []
        for point in corrupted_points:
            for neighbor in point.neighbors:
                if (not neighbor.corrupted) and (not neighbor.neighbor_to_corrupted):
                    points_to_expand.append(neighbor)
                    neighbor.neighbor_to_corrupted = 1
        
        
        t0 = time()
        n_points_to_expend = len(points_to_expand)
        while (n_points_to_expend > 0):
            dn = points_to_expand[0].expand(points_to_expand)
            n_points_to_expend += dn 
            points_to_expand.pop(0)
            n_points_to_expend -= 1
        t1 = time()
        self.djikstra()
        print("mesh : ",t1-t0)
        print("djikstra : ",time()-t1)
        print("---")
        
    def display_djikstra(self):
        # for b in self.beziers:
        #     b.delete()
        # self.beziers = []
        # for i in range(2,len(self.djikstra_path)):
        #     self.beziers.append(Beziers([self.djikstra_path[i],self.djikstra_path[i-1],self.djikstra_path[i-2]],self))
        #     self.beziers[i-2].display()
        for line in self.djikstra_lines:
            self.canvas.delete(line)
        if self.djikstra_path:
            for i in range(1,len(self.djikstra_path)):
                xa,ya = np.array(self.djikstra_path[i-1].pos)*SCREEN_SIZE/10
                xb,yb = np.array(self.djikstra_path[i].pos)*SCREEN_SIZE/10
                ya = SCREEN_SIZE - ya
                yb = SCREEN_SIZE - yb
                self.djikstra_lines.append(self.canvas.create_line(xa, ya, xb, yb,fill = "green" ,width = 3))
            xa,ya = np.array(self.djikstra_path[0].pos)*SCREEN_SIZE/10
            xb,yb = np.array(self.end.pos)*SCREEN_SIZE/10
            ya = SCREEN_SIZE - ya
            yb = SCREEN_SIZE - yb
            self.djikstra_lines.append(self.canvas.create_line(xa, ya, xb, yb,fill = "green", width = 3))
    def display(self):
        t0 = time()
        for line in self.lines:
            self.canvas.delete(line)
        for obstacle in self.obstacles:
            xa,ya = np.array(obstacle.pos_a)*SCREEN_SIZE/10
            xb,yb = np.array(obstacle.pos_b)*SCREEN_SIZE/10
            ya = SCREEN_SIZE - ya
            yb = SCREEN_SIZE - yb
            self.lines.append(self.canvas.create_line(xa, ya, xb, yb,fill = "red"))
        # for point in self.mesh:
        #     x,y = np.array(point.pos)*SCREEN_SIZE/10
        #     y = SCREEN_SIZE - y
        #     for neighbor in point.neighbors:
        #         xn,yn = np.array(neighbor.pos)*SCREEN_SIZE/10
        #         yn = SCREEN_SIZE - yn
        #         self.lines.append(self.canvas.create_line(x, y, xn, yn,fill = "black"))
        self.display_djikstra()
        print("display : ",time()-t0)
    def djikstra(self):
        for point in self.mesh:
            point.djikstra_dist = float("inf")
            point.djikstra_previous = None
            point.djikstra_validated = 0
        over = 0
        points_computed = [self.start]
        self.start.djikstra_dist = 0
        n_points_computed = 1
        while n_points_computed > 0 and not over:
            distances = [point.djikstra_dist for point in points_computed]
            min_idx = distances.index(min(distances))
            chosen_point = points_computed[min_idx]
            if chosen_point.get_distance(self.end)<chosen_point.radius:
                over = 1
            else:
                chosen_point.djikstra_validated = 1
                dist_of_chosen_point = distances[min_idx]
                points_computed.pop(min_idx)
                n_points_computed -= 1
                for neighbor in chosen_point.neighbors:
                    if not neighbor.djikstra_validated:
                        if not neighbor.djikstra_dist<float("inf"):
                            points_computed.append(neighbor)
                            n_points_computed += 1
                        if neighbor.djikstra_dist > dist_of_chosen_point + chosen_point.get_djikstra_distance(neighbor):
                            neighbor.djikstra_dist = dist_of_chosen_point + chosen_point.get_djikstra_distance(neighbor)
                            neighbor.djikstra_previous = chosen_point
        
        self.djikstra_path = []
        if over:
            point = chosen_point
            while point.pos != self.start.pos:
                self.djikstra_path.append(point)
                point = point.djikstra_previous
            self.djikstra_path.append(point)
                     
            
class Beziers():
    def __init__(self,points,world):
        self.points = points
        self.line = None
        self.canvas = world.canvas
    def get_point(self,t):
        sx = 0
        sy = 0
        i = 0
        n = len(self.points)
        for point in self.points:
            x,y = np.array(point.pos)*SCREEN_SIZE/10
            y = SCREEN_SIZE-y
            fact = special.binom(n-1,i)
            sx += fact*x*(t**i)*((1-t)**((n-1)-i))
            sy += fact*y*(t**i)*((1-t)**((n-1)-i))
            i += 1
        return([sx,sy])
    def delete(self):
        if self.line:
             self.canvas.delete(self.line)
    def display(self):
        if self.line:
            self.canvas.delete(self.line)
        xy = []
        for t in np.linspace(0,1,10*len(self.points)):
            nx,ny = self.get_point(t)
            xy.append(nx)
            xy.append(ny)
        self.line = self.canvas.create_line(xy, fill='red')
        
            

obstacles = []
obstacles.append(Obstacle([0,0],[10,0]))
obstacles.append(Obstacle([10,0],[10,10]))
obstacles.append(Obstacle([10,10],[0,10]))
obstacles.append(Obstacle([0,10],[0,0]))

obstacles.append(Obstacle([0,5],[4,5]))
obstacles.append(Obstacle([4,5],[4,4]))
obstacles.append(Obstacle([4,3],[4,0]))

obstacles.append(Obstacle([0,7],[3,7]))
obstacles.append(Obstacle([4,7],[4,10]))

obstacles.append(Obstacle([6,10],[6,9]))
obstacles.append(Obstacle([6,8],[6,4]))
obstacles.append(Obstacle([6,3],[6,0]))
obstacles.append(Obstacle([6,5],[10,5]))

my_world = World(obstacles,[1,1],[9,9])

my_robot = Robot(my_world)
t1 = threading.Thread(target=my_robot.update)
t1.start()
my_world.display()
mainloop()

            