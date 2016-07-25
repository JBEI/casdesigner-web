# -*- mode: ruby -*-
# vi: set ft=ruby :

Vagrant.configure(2) do |config|
  config.vm.box = "debian/jessie64"
  config.vm.network "forwarded_port", guest: 8000, host: 7777

  config.vm.provision "shell", inline: <<-SHELL
    sudo apt-get install -y python3-pip 
    sudo apt-get install -y python3-dev # for Pandas
    sudo apt-get install -y pkg-config libfreetype6-dev libpng12-dev # for matplotlib
    sudo pip3 install jupyter biopython intermine pandas matplotlib xlrd
  SHELL
end
