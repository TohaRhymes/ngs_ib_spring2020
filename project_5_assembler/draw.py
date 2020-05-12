from graphviz import render
from graphviz import Source


def visualize(name):
    render('dot', 'png', name)
    # To render an existing file in a notebook
    Source.from_file(name)


if __name__ == '__main__':
    name_save = 'graphs/graph (3).dot'
    visualize(name_save)
