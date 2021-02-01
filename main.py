import re

from kivy.app import App
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.label import Label
from kivy.uix.textinput import TextInput
from kivy.uix.button import Button
from kivy.graphics import *

'''
The homologous_suffixes dictionary shows the meaning of various endings to compounds by describing them as adjacency 
lists which can then be appended to the main adjacency list. Each list is accompanied by a numerical value that 
signifies the bond strength needed to attach to a carbon. This is useful as groups that require more bonds than
available would be unable to connect.

'ane', 'ene' and 'yne' are unique in that they don't mean more atoms need to be added, but that the bonds on the carbon
atoms present need to be modified. In this prototype, this functionality hasn't been added yet, so their lists remain
empty until a system is devised.
'''

homologous_suffixes = {
    'ene': ({}, 0),
    'yne': ({}, 0),
    'ol': ({0: ('O', [(-1, 1)])}, 1),
    'amine': ({0: ('N', [(-1, 1)])}, 1),
    'al': ({0: ('O', [(-1, 2)])}, 2),
    'one': ({0: ('O', [(-1, 2)])}, 2),
    'oic acid': ({0: ('O', [(-1, 2)]), 1: ('O', [(-1, 1)])}, 3),
    'amide': ({0: ('O', [(-1, 2)]), 1: ('N', [(-1, 1)])}, 3)
}

valency = {
    'C': 4,
    'O': 2,
    'N': 3,
    'Cl': 1,
    'Br': 1
}

regex_hmlg_suff = '|'.join(homologous_suffixes.keys())

numbers = dict(nonacont=90, octacont=80, heptacont=70, hexacont=60, pentacont=50, tetracont=40, triacont=30, cos=20,
               dec=10, non=9, oct=8, hept=7, hex=6, pent=5, but=4, tetr=4, prop=3, tr=3, eth=2, d=2, meth=1, un=1,
               hen=1)

regex_numbers = '|'.join(numbers.keys())

'''
The following two functions have repeated applications throughout the process of producing the adjacency list, so have
been declared.

word_to_num takes the written name for a certain number (e.g. dodeca meaning twelve) and outputs the numerical value.
This is useful when finding the number of carbons in the main chain (e.g. dodecane meaning twelve carbons) or the number
of occurrences ofa sub-group (e.g. dodecol meaning twelve alcohol groups).

bond_strength_total calculates the total strength of the bonds around an atom to see if extra neighbours can be added.
Here, a single bond has strength 1, a double has strength 2, and a triple has strength 3. Each carbon can have a max
bond strength of four. If adding a group would make this above four, the group will not be added. 
'''


def word_to_num(string):
    if not isinstance(string, str):
        raise TypeError('Only strings are allowed')
    string = re.sub(rf'[aio]+$', '', string)
    if string in numbers.keys():
        return numbers[string]
    string_end = re.search(rf'({regex_numbers})$', string)
    if string_end:
        string = re.sub(rf'{string_end.group()}$', '', string)
        return numbers[string_end.group()] + word_to_num(string)
    return 0


def bond_strength_total(atom):
    orbital_sum = 0
    for i in atom[1]:
        orbital_sum += i[1]
    return orbital_sum


def c_search(count, x, y, x_fac, y_fac, molecule, visited):
    x += x_fac * 10
    y += y_fac * 10
    molecule[count].append([x, y])
    visited[count] = [x, y]
    for adj_node in molecule[count][1]:
        if adj_node[0] not in visited.keys():
            while [x + x_fac * 10, y + y_fac * 10] in visited.values():
                x_fac, y_fac = -y_fac, x_fac
            molecule, visited = c_search(adj_node[0], x, y, x_fac, y_fac, molecule, visited)
    return molecule, visited


class OrganicCompound:
    def __init__(self, name):
        self.name = name
        self.structure = {}  # this will become the adjacency matrix that represents the whole compound
        self.split_name = re.split('yl', self.name)  # isolates all of string after final 'yl'
        self.parent_chain = self.split_name.pop()  # takes the central part of the compound to digest first
        self.func_group = re.search(rf'({regex_hmlg_suff})$', self.parent_chain)  # finds series name at end of string
        if self.func_group:
            self.func_group = self.func_group.group()
            self.parent_chain = re.sub(rf'({self.func_group})$', '', self.parent_chain)  # removes series from string
        else:
            self.func_group = 'ane'  # defaults to alkane if not specified
        self.func_freq = re.search(rf'((({regex_numbers})a?i?o?)+)$', self.parent_chain)
        # looks for number of occurrences of series
        if self.func_freq:
            self.parent_chain = re.sub(rf'{self.func_freq.group()}$', '', self.parent_chain)  # removes num from string
            self.func_freq = word_to_num(self.func_freq.group())
        else:
            self.func_freq = 1  # defaults to 1 if not specified
        self.func_num = re.search(rf'-(\d|,)+-$', self.parent_chain)  # looks for positions of groups
        if self.func_num:
            self.parent_chain = re.sub(rf'{self.func_num.group()}$', '', self.parent_chain)  # removes pos from string
            self.func_num = re.split(',', self.func_num.group().replace('-', ''))  # splits pos list
        else:
            self.func_num = [1]
        while len(self.func_num) < self.func_freq:
            self.func_num += [1]
        self.parent_chain = re.sub(rf'a?n?$', '', self.parent_chain)  # removes useless characters often found at end
        self.chain_len = re.search(rf'((({regex_numbers})a?i?o?)+)$', self.parent_chain)  # looks for chain length
        if self.chain_len:
            self.parent_chain = re.sub(rf'{self.chain_len.group()}$', '', self.parent_chain)
            self.chain_len = word_to_num(self.chain_len.group())

    def prod_adj_matrix(self):

        for i in range(self.chain_len):
            self.structure[i] = ['C', []]  # adds carbon atom to adj matrix
            if i:  # if not the first carbon in the chain
                if str(i) in self.func_num:
                    if self.func_group == 'ene':
                        strength = 2
                    elif self.func_group == 'yne':
                        strength = 3
                else:
                    strength = 1
                self.structure[i - 1][1].append([i, strength])  # declares bond from previous carbon to it
                self.structure[i][1].append([i - 1, strength])  # declares bond from it to previous carbon

        group = homologous_suffixes[self.func_group][0]  # retrieves adj list of sub group to append to main adj list
        for j in range(self.func_freq):  # for each occurrence of functional group
            orbital_sum = 5
            while orbital_sum > 4 and int(self.func_num[j]) <= self.chain_len:
                orbital_sum = bond_strength_total(self.structure[int(self.func_num[j]) - 1]) + \
                              homologous_suffixes[self.func_group][1]
                if orbital_sum > 4:
                    self.func_num[j] = int(self.func_num[j]) + 1
            if int(self.func_num[j]) <= self.chain_len:
                structure_len = len(self.structure)  # stores length of adj list before to increment positions
                for k in group:  # for each atoms in the group's adj list
                    bonds = []  # rewrites the bond information
                    for m in group[k][1]:  # for each bond from the atom
                        if m[0] < 0:  # -1 refers to a bond to the main adj list
                            bonds.append(
                                [int(self.func_num[j]) - 1, m[1]])  # substitutes in given position of bond to main
                            self.structure[int(self.func_num[j]) - 1][1].append([structure_len + k, m[1]])
                        else:  # otherwise bond is entirely within the sub group
                            bonds.append([structure_len + m[0], m[1]])
                    self.structure[structure_len + k] = [group[k][0], bonds]
                    '''
        current_len = len(self.structure)
        for atom in range(current_len):
            if self.structure[atom][0] == 'C':
                orbital_sum = 0
                while orbital_sum < 4:
                    orbital_sum = bond_strength_total(self.structure[atom])
                    if orbital_sum < 4:
                        structure_len = len(self.structure)
                        self.structure[structure_len] = ('H', [(atom, 1)])
                        self.structure[atom][1].append((structure_len, 1))
                        '''
        print(self.structure)

    def assign_positions(self):
        # assigns each node a position in a 2D space. 0:['C',[(1,1)]] -> 0:['C',[(1,1),(2,1),(3,1),(4,1)],[10,0]]
        self.structure, visited = c_search(0, 0, 1, 1, 0, self.structure, {})
        print(self.structure)
        temp_molecule = [i for i in self.structure.keys()]
        # loops through every node in self.structure
        print('o')
        for node in temp_molecule:
            free_valency = valency[self.structure[node][0]] - bond_strength_total(self.structure[node])
            # print(self.structure)
            node_x, node_y = self.structure[node][2][0], self.structure[node][2][1]
            x, y = 0, 1
            neighbours = [self.structure[i[0]][2] for i in self.structure[node][1]]

            for i in range(free_valency):
                position = len(self.structure.keys())
                while [node_x + 3.5 * x, node_y - 3.5 * y] in visited.values() or [node_x + 10 * x,
                                                                                   node_y - 10 * y] in neighbours:
                    '''or [node_x + 10 * x, node_y - 10 * y] in visited.values():'''  # something's afoot here
                    x, y = y, -x
                self.structure[node][1].append((position, 1))
                self.structure[position] = ['H', [(node, 1)], [node_x + 3.5 * x, node_y - 3.5 * y]]
                visited[position] = [node_x + 3.5 * x, node_y - 3.5 * y]
        print('k')
        min_x = min([i[2][0] for i in self.structure.values()])
        max_x = max([i[2][0] for i in self.structure.values()])
        min_y = min([i[2][1] for i in self.structure.values()])
        max_y = max([i[2][1] for i in self.structure.values()])
        x_range = (max_x - min_x)
        x_ratio = 640 / x_range
        y_range = (max_y - min_y)
        y_ratio = 410 / y_range
        if x_ratio < y_ratio:
            y_ratio = x_ratio
            min_y -= (410 / y_ratio - y_range) / 2
        else:
            x_ratio = y_ratio
            min_x -= (640 / x_ratio - x_range) / 2
        for node in self.structure:
            self.structure[node][2] = [80 + x_ratio * (self.structure[node][2][0] - min_x),
                                       80 + y_ratio * (self.structure[node][2][1] - min_y)]
        print(self.structure)
        return self.structure


'''
The above section takes the name of a compound and produces an adjacency list that describes its structure. In the
finished program, this will be converted into a graph (with nodes and edges) that acts as the displayed formula for the
compound. However, in the prototype, only the adjacency list is produced.

The below section of the code produces the GUI, which facilitates the input of compound names and will, in the finished
version, display an image of the graph produced. Currently, it can be used to input the name and states 'Look at 
console' when the button is pressed. This is to prove that the window can be updated and the input can be processed
'''


class MainScreen(BoxLayout):

    def __init__(self, **kwargs):
        """
        def comp_conv(instance):  # function triggered by pressing comp_conf button
            compound = OrganicCompound(self.comp_name.text)  # retrieves and processes given compound name
            self.comp_structure.text = 'Look at console'
            compound.prod_adj_matrix()  # prints adjacency list in console
            print(compound.assign_positions())"""

        self.subscreen = BoxLayout(orientation='horizontal', size_hint=(1, .05))  # top row of the GUI
        super(MainScreen, self).__init__(**kwargs)
        self.subscreen.add_widget(Label(text='Compound Name', size_hint=(0.17, 1)))
        self.comp_name = TextInput(multiline=False, size_hint=(0.78, 1))  # input for name of compound
        self.subscreen.add_widget(self.comp_name)
        self.comp_conf = Button(text='Go!', size_hint=(0.05, 1))  # button that triggers comp_conv function
        self.comp_conf.bind(on_press=self.comp_conv)
        self.subscreen.add_widget(self.comp_conf)
        self.orientation = 'vertical'
        self.comp_structure = Label(text='', size_hint=(1, .95))  # blank label that is updated on button press
        self.add_widget(self.subscreen)
        with self.comp_structure.canvas:
            Rectangle(pos=[0, 0], size=[800, 570])
        self.add_widget(self.comp_structure)

    def comp_conv(self, spare):
        if self.comp_name.text == '':
            return
        compound = OrganicCompound(self.comp_name.text)
        compound.prod_adj_matrix()
        layout = compound.assign_positions()
        print(layout)
        rad = max(layout[1][2][0] - layout[0][2][0], layout[1][2][1] - layout[0][2][1]) / 5
        print(rad)
        with self.comp_structure.canvas:
            Rectangle(pos=[0, 0], size=[800, 570])
            for i in layout.keys():
                for j in layout[i][1]:
                    x = 0
                    Color(x, x, x)
                    if j[0] > i:
                        if j[1] > 2:
                            Line(points=layout[i][2] + layout[j[0]][2], width=rad * 1.25)
                            x = abs(x - 1)
                            Color(x, x, x)
                        if j[1] > 1:
                            Line(points=layout[i][2] + layout[j[0]][2], width=rad * 0.75)
                            x = abs(x - 1)
                            Color(x, x, x)
                        Line(points=layout[i][2] + layout[j[0]][2], width=rad / 4)
                Color(0, 0, 0)
                if layout[i][0] == 'H':
                    spare_rad = rad
                    rad = rad * 0.625
                Ellipse(pos=[layout[i][2][0] - rad - rad / 4, layout[i][2][1] - rad - rad / 4],
                        size=[2.5 * rad, 2.5 * rad])
                if layout[i][0] == 'C':
                    Color(0.3, 0.3, 0.3)
                elif layout[i][0] == 'O':
                    Color(0.75, 0, 0)
                elif layout[i][0] == 'H':
                    Color(1, 1, 1)
                    rad -= 0.125 * spare_rad
                elif layout[i][0] == 'N':
                    Color(0.5, 0.5, 1)
                Ellipse(pos=[layout[i][2][0] - rad, layout[i][2][1] - rad], size=[2 * rad, 2 * rad])
                if layout[i][0] == 'H':
                    rad = spare_rad


class MyApp(App):

    def build(self):
        self.borderless = True
        self.title = 'Compound Sketcher'
        return MainScreen()


if __name__ == '__main__':
    MyApp().run()

'''Examples to test the program:
    propan-2-ol
       H  H  H
       |  |  |
    H--C--C--C--H
       |  |  |
       H  O  H
          H
    propan-2,2,2-triol
    INVALID AS ONLY TWO alcohol GROUPS CAN FIT ON CENTRAL CARBON. WILL DEFAULT TO
    propan-2,2,3-triol [SAME AS propan-1,2,2-triol]
          H
       H  O  H
       |  |  |
    H--C--C--C--OH
       |  |  |
       H  O  H
          H
    propan-1,3-dioic acid
        H  H  H
        |  |  |
     O==C--C--C==O
        |  |  |
        O  H  O
        H     H



'''
