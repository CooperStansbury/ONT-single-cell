import pandas as pd
import os
import sys


def getFastq(path):
    """A function to import basecalls """
    imageDf = pd.read_csv(path,  comment="#")
    return imageDf

