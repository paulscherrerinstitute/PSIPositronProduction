# How to build VSCode from source on PSI RHEL7

I am a big fan of the Microsoft Visual Studio Code editor.
At today (11.08.2021), no official rpm packages for Red Hat Enterprise Linux 7.9 (RHEL7) are provided anymore.
In this page we provide specific, step-by-step instructions to build the Microsoft Visual Studio Code editor (VSCode) from source on the mentioned PSI official Linux system.


## Reference page

The reference for the development of VSCode is:

https://github.com/microsoft/vscode

Instructions to build and run from source are at:

https://github.com/microsoft/vscode/wiki/How-to-Contribute


## Specific prerequisites for RHEL7

### Recent GCC toolchain

See [How to build Geant4 on PSI RHEL7](Geant4OnPsiRhel7.md#recent-gcc-toolchain)


### Node.JS, newest versions (e.g. 14.x)

The newest versions are not available in the RHEL7 official repositories.

Add the Node.JS repository with:

```shell
$ curl -sL https://rpm.nodesource.com/setup_14.x | sudo bash -
```

Install Node.JS with:

```shell
$ sudo yum install -y nodejs
```

Verify the installed verison:

```shell
$ node -v
$ v14.17.4
$ npm -v
$ 6.14.14
```

For more details, see
https://computingforgeeks.com/install-node-js-14-on-centos-rhel-linux
.


### Packages fakeroot and rpmdevtools

If you want to create an rpm package (not necessary for a singel user):

```shell
$ sudo yum install fakeroot rpmdevtools
```

For more details, see
https://access.redhat.com/documentation/en-us/red_hat_enterprise_linux/7/html-single/rpm_packaging_guide
.


## Install VSCode

### Build

This sequence of commands is known to work (do not forget to activate devtoolset!):

```shell
$ git clone https://github.com/microsoft/vscode.git
$ cd vscode
$ npm install yarn
$ ./node_modules/yarn/bin/yarn
$ ./node_modules/yarn/bin/yarn run gulp vscode-linux-x64
```

After the last command succesfully finishes, you will find the executable in `../VSCode-linux-x64/bin/code-oss`.

Before being able to run it, it is important to set the correct owner and permissions to `chrome-sandbox`:

```shell
$ cd ../VSCode-linux-x64
$ sudo chown root:root chrome-sandbox
$ sudo chmod 4755 chrome-sandbox
```


### Generate and install rpm

Build the rpm package with:

```shell
$ ./node_modules/yarn/bin/yarn run gulp vscode-linux-x64-build-rpm
```

Install the generated package with:

```shell
$ sudo yum localinstall ./.build/linux/rpm/x86_64/code-oss-*.el7.x86_64.rpm
```


### Start VSCode

```shell
$ code-oss
```

or, for more information:

```shell
$ code-oss --verbose
```

