#!/usr/bin/env python3
# encoding: utf-8
# 2023 Harald Klimach <harald.klimach@dlr.de>

APPNAME = 'treelm'

top = '.'
out = 'build'

def options(opt):
    '''Building options provided by the TreElM-Library.
       Remember, all options can be displayed with waf --help.'''
    opt.recurse('bin')
    opt.load('coco')
    opt.recurse('aotus')
    opt.recurse('tem')


def configure(conf):
    '''Project configuration'''
    import os
    conf.recurse('aotus', 'subconf')
    conf.recurse('bin', 'preconfigure')
    # Initialize the coco preprocessing tool
    if not conf.options.coco_set:
        conf.options.coco_set = 'default.coco'
    conf.load('coco')
    conf.recurse('tem')
    conf.recurse('bin', 'postconfigure')


def build(bld):
    '''Build the treelm project'''
    from revision_module import fill_revision_string
    bld.recurse('bin')
    if not (bld.cmd == 'docu' and bld.env.fordonline):
        bld.recurse('aotus')
    fill_revision_string(bld, subdir='tem')
    bld.recurse('tem')
