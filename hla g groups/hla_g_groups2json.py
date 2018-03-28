import xml.etree.ElementTree
from zipfile import ZipFile
import requests
import io
import json

XML_URL = 'https://github.com/ANHIG/IMGTHLA/raw/Latest/xml/hla_ambigs.xml.zip'

def download_ambigs_xml(url):
    '''
    download the ambiguity zipped xml from github into memory
    '''
    #download
    dl = requests.get(url)
    #convert dl.content to Byte stream
    #convert as ZipFile object
    #read the wanted file
    #return a file like object
    ZipFile(io.BytesIO(dl.content)).extractall()

    
def get_ggroup(ambiguities_xml_file, 
               reverse=False, 
               only_first_allele=False):
    
    #read the given xml file downloaded from imgt
    #geneList   0 = releaseVersion
    hla = xml.etree.ElementTree.parse(ambiguities_xml_file).getroot().getchildren()[1] 
    
    gGroups = {}
    
    #iterate through the file to get all the locus
    for loci in hla.getchildren():
        locus_name = loci.attrib['name']
        #get the locus gGroup = first child element
        locus_gGroup = loci.getchildren()[0]
        
        for gGroup in locus_gGroup.getchildren():
            gGroupName = gGroup.attrib['name'].replace('HLA-', '')
            
            if only_first_allele:
                firstAllele = gGroup.getchildren()[0].attrib['name'].replace('HLA-', '')
            
                if reverse:
                    gGroups[firstAllele] = gGroupName
                else:
                    gGroups[gGroupName] = firstAllele                    
            else:
                alleles = [allele.attrib['name'].replace('HLA-', '') \
                           for allele in gGroup.getchildren()]
                gGroups[gGroupName] = alleles
        
    return gGroups

def write_to_json(json_file):
    '''
    make a json file for hla_ambigs.xml only keeping the G groups
    '''
    with open('hla_ambigs.json', 'w') as fout:
        fout.write(json.dumps(get_ggroup('hla_ambigs.xml')))

        
        
if __name__ == '__main__':
    #download
    print('\n\tDownloading the latest HLA amibiguity list xml file from Github...')
    download_ambigs_xml(XML_URL)
    print('\n\tConverting the list to json...')
    #write to a local json file
    write_to_json('hla_ambigs.json')
    print('\n\thla_ambig.json has been updated to the latest version...\n')