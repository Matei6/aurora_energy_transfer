# -*- coding: utf-8 -*-
"""
Created on Sun Mar 23 18:26:54 2025

@author: matei
"""
import random as rd


def main():
    import turtle

    turtle.clearscreen()
    turtle.up()
    turtle.goto(-330, 0)
    turtle.down()
    turtle.speed(7)
    angle = [0, rd.randint(-55, 55)]

    for _ in range(100):
        turtle.forward(rd.randint(5, 15))
        turtle.left(angle[1] - angle[0])
        angle = [angle[1], rd.randint(-55, 55)]

    turtle.mainloop()
    turtle.done()


if __name__ == "__main__":
    main()