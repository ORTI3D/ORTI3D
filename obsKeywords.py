# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 09:42:04 2014

@author: olive
"""
class Obs:
    def __init__(self):
        self.grpList =['OBS']
        self.groups = {'OBS':['obs.1']}
        self.lines = {
        'obs.1':{'comm':'observation','kw':['OBS'],'type':['arrint'],'cond':''}
        }
