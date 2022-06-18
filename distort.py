import click
import csv as csvlib
import math

import pewpew.laser_events as l
from pewpew import MachineConnection

def calibration_grid(border, points):

    delta = (1.0 - 2 * border) / (points - 1)
    row = list(border + delta * i for i in range(points))
    middle = points // 2
    
    for i,x in enumerate(row):
        for j, y in enumerate(row):
            circle = (j == 0 and i == 0) or (j == middle and (i == 0 or i == points - 1)) 
            
            yield x, y, circle


def generate_toolpath(grid, speed):
    # Figure out the distance between some adjacent points and use that to determine the
    # radius of circles
    (a,b,_),(c,d,_) = grid[0:2]
    radius = 0.25 * math.sqrt((a-c)**2 + (b - d)**2)
    length = 0.25 * radius

    def gen():
        for x,y,circle in grid:
            yield l.line((x - length, y),(x + length, y), speed = speed)
            yield l.line((x, y - length),(x, y + length), speed = speed)
            
            if circle is True:
                yield l.circle((x,y), radius, 0, 2 * math.pi, speed = speed)
                
    return list(l.adjust_delays(gen(), 10.0))

@click.command()
@click.argument("device")
@click.option('--preview/--engrave', default = False, help = "Execute the pattern with the guide laser")
@click.option('--power', default = 0.25, type = click.FloatRange(0.0,1.0), help = "Laser power (0-1.0)")
@click.option('--speed', default = 4, type = click.FloatRange(0.0,20.0), help = "Laser speed for marking circles")
@click.option('--points', default = 11, type = int, help = "Grid size for calibation - must be odd, greater than 1")
@click.option('--border', default = 0.05, type = click.FloatRange(0,0.5), help = "Border width")
@click.option('--csv', default = 'calibration_grid.csv', type = click.File('w'), help = "Output file for calibration grid points")

def cal(device, preview, power, speed, points, border, csv):

    if points < 3 or 0 == points % 2:
        print("Invalid grid size - must be odd, greater than 1")
    
    grid = list(calibration_grid(border, points))
    csvlib.writer(csv).writerows([['X', 'Y', 'Circled']] + grid)
    toolpath = list(generate_toolpath(grid, speed))


    print("Connecting to laser...")
    with MachineConnection(device) as m:

        m.buffered_messages([l.arm(),l.enable_guide()])
        m.wait_until_idle()

        if preview:
            input("Press enter to begin grid preview ")
        else:
            print("Caution! Laser is armed and ready to fire!")
            input("Press enter to begin ")
            m.buffered_messages([l.disable_guide(), l.power(power)])

        m.buffered_messages(toolpath)
        m.buffered_messages([l.disarm()])
        m.wait_until_idle()

        
    
    
if __name__ == '__main__':
    cal()

