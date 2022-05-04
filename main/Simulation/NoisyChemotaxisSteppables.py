import numpy as np

from cc3d.core.PySteppables import *
#change volume
#change amount of noise in env

#global variable selection

def distance(coord_1, coord_2, p=None):

    if p:
        print(coord_1 - coord_2)
    x = np.linalg.norm(coord_1 - coord_2)
    return x


class NoisyChemotaxisSteppable(SteppableBasePy):
    def __init__(self, frequency=1):

        SteppableBasePy.__init__(self, frequency)

    def start(self):
        self.build_wall(self.WALL)
        # any code in the start function runs before MCS=0
        # Set Parameter values here and create persistent variables

        
        self.plot_win = self.add_new_plot_window(title='Cell Trajectory',
                                                 x_axis_title='X-Location',
                                                 y_axis_title='Y-Location', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False)
        
        self.plot_win.add_plot("CellLoc", style='Lines', color='red', size=5)

        
        self.plot_win_2 = self.add_new_plot_window(title='Selected Parameter vs Speed',
                                                 x_axis_title='Parameter Value',
                                                 y_axis_title='Speed / MCS', x_scale_type='linear', y_scale_type='linear',
                                                   grid=False, config_options={"legend":True})
        
        self.plot_win_2.add_plot("X_Speed", style='Dots', color='Red', size=10)
        self.plot_win_2.add_plot("Y_Speed", style='Dots', color='Blue', size=10)
        self.plot_win_2.add_plot("Tot_Speed", style='Dots', color='Purple', size=10)

        
        


        test = False
        self.shared_steppable_vars["prob"]= 0.01 ## Probability of there being holes in the gradient
        self.shared_steppable_vars["runtime"]= 2000 ## Probability of there being holes in the gradient
        self.field = CompuCell.getConcentrationField(self.simulator, "glucose")

        self.noise = 1.2 ## How much the current gradient is divided by [1,].
        at_gradient = 0
        flux = 1 ## Should initialize at 0 once we have a steering panel

        ## Setting the size of the lattice for the field
        self.shared_steppable_vars["lattice_dimensions"] = [120, 60]
        self.x_lattice_size = self.shared_steppable_vars["lattice_dimensions"][0]
        self.y_lattice_size = self.shared_steppable_vars["lattice_dimensions"][1]


        self.shared_steppable_vars["lambda_chemo"] = 800.0
        self.shared_steppable_vars["targetVolume"] = 50.0
        self.shared_steppable_vars["selected_parameter"] = None
        self.shared_steppable_vars["last_measure"] = 0


        ##############################################################################################
        ##############################################################################################
        ##############################################################################################
        ##############################################################################################

        for cell in self.cellList:

            if cell.type == 1:

                cell.dict['X_Speeds'] = []
                cell.dict['Y_Speeds'] = []

                cd = self.chemotaxisPlugin.addChemotaxisData(
                    cell, "glucose"
                )  # Cell and field it responds to

                cd.setLambda(
                    self.shared_steppable_vars["lambda_chemo"]
                )  # Strength of attraction

                # cd.assignChemotactTowardsVectorTypes([self.MEDIUM, self.BACTERIUM])
                cell.lambdaVolume = 5.0  # In this simulation I am scanning lambda volume between min_value and max_value
                cell.targetVolume = self.shared_steppable_vars["targetVolume"]
                r = np.sqrt(cell.targetVolume/np.pi)
                cell.targetSurface= 2*np.pi*r
                cell.lambdaSurface= 0.5


        # patchy gradient
        secrConst = 1  # starting constant

        self.shared_steppable_vars["starting_location"] = None

        self.initialize_season(
            lattice_size=self.shared_steppable_vars["lattice_dimensions"], set_prob= self.shared_steppable_vars["prob"]
        )

    def step(self, mcs):
        # # any code in the start function runs before MCS=0
        # # Set Parameter values here and create persistent variables
        # arguments are (name of the data series, x, y)
        # iterating over all cells in simulation

        self.shared_steppable_vars['current_mcs'] = mcs

        while self.get_steering_param('variable selection')=="Select Variable to Manipulate":
            pass


        for cell in self.cellList:
            if cell.type == 1:

                if mcs > 1:
                    if self.get_steering_param('variable selection')=='volume' and (mcs-self.shared_steppable_vars['last_measure'])<20:
                        continue
                    else:
                        cell.dict['X_Speeds'].append(abs(cell.xCOM - cell.xCOMPrev))
                        cell.dict['Y_Speeds'].append(abs(cell.yCOM - cell.yCOMPrev))
                    

                    if self.get_steering_param('variable selection')=="noise":
                        if (mcs - self.shared_steppable_vars['last_measure']) == self.shared_steppable_vars['runtime']:
                            self.shared_steppable_vars['last_measure'] = mcs
                            # arguments are (name of the data series, x, y)
                            self.plot_win_2.add_data_point("X_Speed", self.shared_steppable_vars['prob'], np.mean(cell.dict['X_Speeds']))
                            self.plot_win_2.add_data_point("Y_Speed", self.shared_steppable_vars['prob'], np.mean(cell.dict['Y_Speeds']))
                            self.plot_win_2.add_data_point("Tot_Speed", self.shared_steppable_vars['prob'], (np.mean(cell.dict['Y_Speeds'])+np.mean(cell.dict['X_Speeds'])))
                            ## Refresh the speeds so that they do not impact the new calculations
                            cell.dict['X_Speeds']=[]
                            cell.dict['Y_Speeds']=[]
                            if self.shared_steppable_vars['prob'] < 1.0:
                                self.shared_steppable_vars['prob']+=0.01

                            else:
                                self.shared_steppable_vars['prob']=0.01

                            self.initialize_season(
                                lattice_size=self.shared_steppable_vars["lattice_dimensions"], set_prob= self.shared_steppable_vars['prob'], static=True
                            )
    
    
                    if self.get_steering_param('variable selection')=="volume":

                        if (mcs - self.shared_steppable_vars['last_measure']-20) == self.shared_steppable_vars['runtime']:
                            self.shared_steppable_vars['last_measure'] = mcs
                            # arguments are (name of the data series, x, y)
                            self.plot_win_2.add_data_point("X_Speed",cell.targetVolume, np.mean(cell.dict['X_Speeds']))
                            self.plot_win_2.add_data_point("Y_Speed", cell.targetVolume, np.mean(cell.dict['Y_Speeds']))
                            self.plot_win_2.add_data_point("Tot_Speed", cell.targetVolume, (np.mean(cell.dict['Y_Speeds'])+np.mean(cell.dict['X_Speeds'])))
                            ## Refresh the speeds so that they do not impact the new calculations
                            cell.dict['X_Speeds']=[]
                            cell.dict['Y_Speeds']=[]
                            if cell.targetVolume < 300.0:
                                cell.targetVolume+=25
                            else:
                                cell.targetVolume = 25


                    

                if not mcs%10:
                    self.plot_win.add_data_point("CellLoc", cell.xCOM, cell.yCOM)

                if (cell.xCOM < self.shared_steppable_vars["starting_location"][0] + 5) and (cell.xCOM > self.shared_steppable_vars["starting_location"][0] - 5):
                    self.initialize_season(lattice_size=self.shared_steppable_vars["lattice_dimensions"], set_prob=self.shared_steppable_vars['prob'])
                    self.plot_win.erase_all_data()
                    



    def finish(self):
        """
        Finish Function is called after the last MCS
        """

    def on_stop(self):
        # this gets called each time user stops simulation
        return

    def add_steering_panel(self):

        # This should allow them to choose what they want to change.
        self.add_steering_param(name='initial volume', val=25, enum=[(x+1)*25 for x in range(12)],
                                    widget_name='combobox')
        self.add_steering_param(name='initial noise', val=0.00, enum=[(x)*0.01 for x in range(101)],
                                    widget_name='combobox')
        self.add_steering_param(name='variable selection', val='Select Variable to Manipulate', enum=['volume','noise'],
                                    widget_name='combobox')

        
        # self.add_steering_param(name='target_vol', val=25, min_val=25, max_val=300, widget_name='slider')
        # self.add_steering_param(name='percent_noise', val=0.01, min_val=0.0, max_val=0.45, widget_name='slider')

    def initialize_season(self, slope=0.01, baseline=1, lattice_size=[120, 60], set_prob=0.1, static=False):

        # if first starting locations, initialize empty list and pick random location
        secrConst = slope
        if self.shared_steppable_vars["starting_location"] is None:
            self.shared_steppable_vars["starting_location"] = np.zeros(
                shape=(2)
            )  # Initialize empty coordinates
            self.shared_steppable_vars["starting_location"][0] = 110  # Setting the x-coordinate
            self.shared_steppable_vars["starting_location"][1] = 30  # Setting the y-coordinate
        else:
            if not static:
                self.shared_steppable_vars["starting_location"][0] = abs(self.shared_steppable_vars["starting_location"][0] - 120)  # Setting the x-coordinate

        bL_corner = np.array([0, 0])
        uR_corner = np.array([lattice_size[0], lattice_size[1]])
        uL_corner = np.array([0, lattice_size[1]])
        bR_corner = np.array([lattice_size[0], 0])
        distances = [
            distance(self.shared_steppable_vars["starting_location"], bL_corner),
            distance(self.shared_steppable_vars["starting_location"], uR_corner),
            distance(self.shared_steppable_vars["starting_location"], uL_corner),
            distance(self.shared_steppable_vars["starting_location"], bR_corner),
        ]
        self.shared_steppable_vars["max_distance"] = max(distances)
        max_distance = self.shared_steppable_vars["max_distance"]

        # iterate through all pixels, pixel is either secrConstant*distance or zero
        for x, y, z in self.everyPixel(1, 1, 1):
            tracker = self.cellField[x, y, z]
            d = distance(
                np.array([x, y]), self.shared_steppable_vars["starting_location"]
            )
            dist = abs(max_distance - d)
            outcome = np.random.binomial(1, 1 - set_prob)
            if tracker:
                if outcome:
                    self.field[x, y, z] = dist * slope + baseline
                else:
                    #self.field[x, y, z] = 0
                    self.field[x, y, z] = (dist * slope + baseline)/self.noise
            else:
                if outcome:
                    self.field[x, y, z] = dist * slope + baseline
                else:
                    #self.field[x, y, z] = 0
                    self.field[x, y, z] = (dist * slope + baseline)/self.noise

    def process_steering_panel_data(self):
        #global variable selection

        if self.shared_steppable_vars['current_mcs'] == 0:
            selection = self.get_steering_param('variable selection')
            tot_vol = self.get_steering_param('initial volume')
            self.shared_steppable_vars['targetVolume'] = tot_vol
            percent_noise= self.get_steering_param('initial noise')

        if self.get_steering_param('variable selection') == 'volume':
            if (self.shared_steppable_vars['current_mcs'] - self.shared_steppable_vars['last_measure'] > self.shared_steppable_vars['runtime']) and target_vol < 300.0:
                tot_vol += 25.0
                self.shared_steppable_vars['targetVolume'] = tot_vol
        else:
            if (self.shared_steppable_vars['current_mcs'] - self.shared_steppable_vars['last_measure'] > self.shared_steppable_vars['runtime']) and percent_noise < 1.0:
                percent_noise += 0.01
                self.shared_steppable_vars['prob'] = percent_noise



        for cell in self.cell_list:

            cell.targetVolume = tot_vol
            r = np.sqrt(cell.targetVolume/np.pi)
            cell.targetSurface = 2*np.pi*r


        if percent_noise!= self.shared_steppable_vars['prob']:
            self.shared_steppable_vars['prob'] = percent_noise
            self.shared_steppable_vars['last_measure'] = self.shared_steppable_vars['current_mcs']
            self.initialize_season(
                lattice_size=self.shared_steppable_vars["lattice_dimensions"], set_prob= percent_noise, static=True
            )

