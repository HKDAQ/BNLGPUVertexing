
paramfile="parameters.txt"

cd ~/files/
#list_of_runs=`ls -d eminus*14per`
list_of_runs=`ls -d eminus*40per`
cd -


for file in ${list_of_runs}; do

    echo RUNNING FILE ${file}

    E=`echo ${file} | awk -F "_" '{print $2}' | awk -F "MeV" '{print $1}' `
    coverage=`echo ${file} | awk -F "_" '{print $6}' | awk -F "per" '{print $1}' `

    ln -s ~/files/${file}/*txt .


### thresholds for 14%
#    thresholds='50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100'

### thresholds for 40%
#    thresholds='150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 250 251 252 253 254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270 271 272 273 274 275 276 277 278 279 280 281 282 283 284 285 286 287 288 289 290 291 292 293 294 295 296 297 298 299 300'


    ./daq_nhits
	
    mv all_hits_emerald*.txt ~/files/${file}/
	
    rm all_hits_*txt
    rm all_pmts.txt
    rm detector.txt

done
