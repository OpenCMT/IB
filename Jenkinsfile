pipeline {

    agent {
        docker {
            image 'infnpd/cmt-environment:latest-centos7'
            args '-u 0:0'
            label DOCKER
        }
    }

    parameters {
        string(name: 'repourl', defaultValue: 'http://artifacts.pd.infn.it/packages/MUOTOM/rpms/repos/centos7/cmt.repo', description: 'Repository URL')
    }

    stages {

        stage('Build') {

            steps {
                sh "wget -O /etc/yum.repos.d/cmt.repo ${params.repourl}"
                sh "yum -y install cmt-ulib-devel"
                sh "mkdir build && cd build && cmake -D ULIB_USE_QT5:BOOL=OFF .. && make"
            }

        }

    }

}
