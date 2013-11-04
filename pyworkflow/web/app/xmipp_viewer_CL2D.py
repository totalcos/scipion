import os
from views_util import *

############## VIEWER XMIPP CL2D ##############
def viewerCL2D(request, protocol, protocolViewer):
    ioDict = {}
    
#    if protocolViewer.doShowLastLevel:
    typeUrl, url = viewLevelFiles(request, protocol, protocolViewer)
    ioDict[typeUrl]= url
    
    if protocolViewer.doShowClassHierarchy:
        typeUrl, url = viewClassHierarchy(request, protocol, protocolViewer)
        ioDict[typeUrl]= url
    
    return ioDict

def viewLevelFiles(request, protocol, protocolViewer):
    fnSubset = protocolViewer._getSubset()
    levelFiles = protocol._getLevelMdFiles(fnSubset)
    
    if levelFiles:
        levelFiles.sort()
        lastLevelFile = levelFiles[-1]
        if protocolViewer.doShowLastLevel:
            return "showj", "/visualize_object/?path="+lastLevelFile
        else:
            if protocolViewer.showSeveralLevels.empty():
                return 'error','Please select the levels that you want to visualize.'
            else:
                listOfLevels = []
                try:
                    listOfLevels = protocolViewer._getListFromRangeString(protocolViewer.showSeveralLevels.get())
                except Exception:
                    return 'error','Invalid levels range.'
                    
                files = "";
                for level in listOfLevels:
                    fn = protocol._getExtraPath("level_%02d/level_classes%s.xmd"%(level,fnSubset))
                    if os.path.exists(fn):
                        files += "/visualize_object/?path="+ fn + "-"
                    else:
                        return 'error','Level %s does not exist.' % level
                if files != "":
                    files = files.split("-")
                    files.pop()
                    return "showjs", files
                    

def viewClassHierarchy(request, protocol, protocolViewer):
    fnSubset = protocolViewer._getSubset()
    fnHierarchy = protocol._getExtraPath("classes%s_hierarchy.txt" % fnSubset)
    if os.path.exists(fnHierarchy):
        html = textfileViewer(fnHierarchy, [fnHierarchy])
    return "html", html
