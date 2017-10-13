RELEASE=`lsb_release -sc`


# Add qgis repository
echo -ne " Adding the Ubuntu GIS unstable repository ..." &&
add-apt-repository -y ppa:ubuntugis/ubuntugis-unstable

if grep -q "qgis.org/ubuntugis" /etc/apt/sources.list;then
	echo " Yeah, you are already a QGIS user, nice!"
else

	echo " "
	echo " Adding the QGIS repository"
	echo " Note: QGIS will not be installed, in order to do so type: "
	echo " sudo apt-get install --yes qgis libqgis-dev"
	echo " "
	echo "deb http://qgis.org/ubuntugis ${RELEASE} main" >> /etc/apt/sources.list
	echo "deb-src http://qgis.org/ubuntugis ${RELEASE} main" >> /etc/apt/sources.list


	wget -O - http://qgis.org/downloads/qgis-2016.gpg.key | gpg --import
	gpg --export --armor 073D307A618E5811 | apt-key add - 
fi

# Update repository
apt-get update -y


# 3 install packages
#------------------------------------------------------------------

SECONDS=0
echo -ne " Installing GIS/Remote sensing packages ..."
apt-get install --yes --allow-unauthenticated \
											gdal-bin \
											libgdal-dev \
											saga \
											libsaga-dev \
											geotiff-bin \
											libgeotiff-dev \
											spatialite-bin \
											spatialite-gui \
											python-scipy



# Install OTB 6.0
otb=OTB-contrib-6.0.0-Linux64
wget https://www.orfeo-toolbox.org/packages/$otb.run
chmod +x $otb.run
mv $otb.run /usr/local/lib
cd /usr/local/lib
./$otb.run
rm $otb.run
ln -s $otb orfeo
chmod o+rx orfeo/*.sh
chmod o+rx orfeo/otbenv.profile
cd -
echo "PYTHONPATH=/usr/local/lib/OTB-contrib-6.0.0-Linux64/lib/python" >> /etc/environment
echo "PATH=${PATH}:/usr/local/lib/orfeo/bin" >> /etc/environment

# Install python RIOS library

wget https://bitbucket.org/chchrsc/rios/downloads/rios-1.4.4.tar.gz
tar -xzf rios-1.4.4.tar.gz
mkdir RIOSInstall
cd rios-1.4.4
python setup.py install --prefix=../RIOSInstall
mv ../RIOSInstall /usr/local/lib
echo "PYTHONPATH=/usr/local/lib/RIOSInstall/lib/python2.7/site-packages" >> /etc/environment
echo "PATH=${PATH}:/usr/local/lib/RIOSInstall/bin" >> /etc/environment