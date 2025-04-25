import json

class ObsBlock:
    """
    A class to represent an observation block.
    """

    def __init__(self, target_name, coversheet_id):
        """
        Initialize the ObsBlock with an ID and name.

        :param obs_block_id: The ID of the observation block.
        :param obs_block_name: The name of the observation block.
        """
        self.target_name = target_name
        self.coversheet_id = coversheet_id
        self.ra = None
        self.dec = None
        self.pm_ra = None
        self.pm_dec = None
        self.pri = None
        self.instrument = None
        self.exposure_time = 0.0
        self.number_of_exposures = 0
        
    def __repr__(self):
        return f"ObsBlock(name={self.target_name}, coversheet_id={self.coversheet_id})"
    
    def to_dict(self):
        """
        Convert the ObsBlock instance to a dictionary.

        :return: A dictionary representation of the ObsBlock.
        """
        return {
            "target_name": self.target_name,
            "coversheet_id": self.coversheet_id,
            "ra": self.ra,
            "dec": self.dec,
            "pm_ra": self.pm_ra,
            "pm_dec": self.pm_dec,
            "pri": self.pri,
            "instrument": self.instrument,
            "exposure_time": self.exposure_time,
            "number_of_exposures": self.number_of_exposures
        }
    def from_dict(self, data):
        """
        Populate the ObsBlock instance from a dictionary.
        :param data: A dictionary containing the ObsBlock data.
        """ 
        self.target_name = data.get("target_name")
        self.coversheet_id = data.get("coversheet_id")
        self.ra = data.get("ra")
        self.dec = data.get("dec")
        self.pm_ra = data.get("pm_ra")
        self.pm_dec = data.get("pm_dec")
        self.pri = data.get("pri")
        self.instrument = data.get("instrument")
        self.exposure_time = data.get("exposure_time")
        self.number_of_exposures = data.get("number_of_exposures")

    def to_json(self):
        """
        Convert the ObsBlock instance to a JSON string.

        :return: A JSON string representation of the ObsBlock.
        """
        return json.dumps(self.to_dict(), indent=4)
    
    def from_json(self, json_str):
        """
        Populate the ObsBlock instance from a JSON string.

        :param json_str: A JSON string containing the ObsBlock data.
        """
        data = json.loads(json_str)
        self.from_dict(data)

    def save_to_file(self, filename):
        """
        Save the ObsBlock instance to a file in JSON format.

        :param filename: The name of the file to save the ObsBlock to.
        """
        with open(filename, 'w') as f:
            f.write(self.to_json())
    def load_from_file(self, filename):
        """
        Load the ObsBlock instance from a file in JSON format.

        :param filename: The name of the file to load the ObsBlock from.
        """
        with open(filename, 'r') as f:
            json_str = f.read()
            self.from_json(json_str)


class hires_obsblock(ObsBlock):
    """
    A class to represent a high-resolution spectrometer observation block.
    """

    def __init__(self, target_name, coversheet_id):
        """
        Initialize the hires_obsblock with an ID and name.

        :param obs_block_id: The ID of the observation block.
        :param obs_block_name: The name of the observation block.
        """
        super().__init__(target_name, coversheet_id)
        self.instrument = "Levy"
        self.iodine = None
        self.decker = None
        self.binning = None
        self.vmag = None
        self.bmv = None
        self.expcount = None
        self.do = None


    def __repr__(self):
        return f"hires_obsblock(name={self.target_name}, coversheet_id={self.coversheet_id})"
    
    def make_scriptobs_line(self, foc, uth, utm, comment=""):
        """
        Create a script observation line for the high-resolution spectrometer.
        """
        line = f"{self.target_name} {self.ra} {self.dec} pm_ra={self.pm_ra} pm_dec={self.pm_dec} "
        line += f"I2={self.iodine} lamp=none uth{uth} utm={utm} "
        line += f"vmag={self.vmag} texp={self.exposure_time} nexp={self.number_of_exposures} foc={foc} "
        line += f"decker={self.decker} binning={self.binning} "
        line += f"do={self.do} owner={self.coversheet_id} "
        line += f"{comment}"

        return line