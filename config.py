import os
import platform
import seq_search
if platform.system() == 'Windows':
    from winreg import *


class Config():
    @staticmethod
    def get_modeller_path():
        try:
            if platform.system() == 'Windows':
                pyprotmodel_key = OpenKey(HKEY_CURRENT_USER, r'SOFTWARE\PyProtModel', 0, KEY_READ)
                [pathVal, regtype] = (QueryValueEx(pyprotmodel_key, 'MODELLER_PATH'))
                CloseKey(pyprotmodel_key)
                return pathVal
            elif platform.system() == 'Linux':
                return os.environ['MODELLER_PATH']
        except:
            seq_search.SeqSearch.modellerConLable.setStyleSheet('color: red')
            seq_search.SeqSearch.modellerConLable.setText('Not connected')

    @staticmethod
    def set_modeller_path(p, v):
        try:
            if platform.system() == 'Windows':
                keyVal = r"SOFTWARE\PyProtModel"
                # create a key in reg if not exist
                if not os.path.exists(keyVal):
                    CreateKey(HKEY_CURRENT_USER, keyVal)
                # set modeller path
                RegistryKey = OpenKey(HKEY_CURRENT_USER, r"SOFTWARE\PyProtModel", 0, KEY_WRITE)
                SetValueEx(RegistryKey, "MODELLER_PATH", 0, REG_SZ, p)
                # set modeller version
                SetValueEx(RegistryKey, "MODELLER_VERSION", 0, REG_SZ, v)
                CloseKey(RegistryKey)
            elif platform.system() == 'Linux':
                os.system("echo 'export MODELLER_PATH=\"%s\"' >> ~/.bashrc" % p)
                os.system("echo 'export MODELLER_VERSION=\"%s\"' >> ~/.bashrc" % v)
                # Popen("source ~/.bashrc", shell=True, executable="/usr/bin/bash")
        except Exception as er:
            print(er)