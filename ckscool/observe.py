"""
Module to plan CKS-Cool observations.

- Calculate needed observing time
- Calculate total open shutter time for observations (anything observed before 2018) doesn't count

- Generate HIRES scripts
- Generate NIRC2 scripts

"""


class ObserveNIRC2(object):
    def generate_nirc2_observing_list():
        pass 

    def expected_observing_time():
        """
        Return amount of on sky NIRC2 LGS time we would need.
        """
        pass 

    def remaining_observing_time():
        """
        Return amount of on sky NIRC2 LGS time we would need.
        """
        pass 

class ObserveHIRES(object):
    def generate_hires_observing_list():
        pass

    def expected_observing_time():
        """
        Return amount of on sky NIRC2 LGS time we would need.
        """
        pass
