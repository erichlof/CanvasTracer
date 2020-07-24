// first, roll my own 3D Vector/Surface/RayTracing Math utility functions ;-)

let a, b, c;
let discrim, rootDiscrim, Q;
let t0 = 0;
let t1 = 0;
let inverseMag = 0;
let ratioIoR = 0;

function vec3(x, y, z)
{
        this.x = x || 0.0;
        this.y = y || 0.0;
        this.z = z || 0.0;     
}

let tempVec = new vec3();
let sN = new vec3();
let L = new vec3();

vec3.prototype.set = function(x, y, z)
{
        this.x = x;
        this.y = y;
        this.z = z;
};

vec3.prototype.copy = function(otherVec)
{
        this.x = otherVec.x;
        this.y = otherVec.y;
        this.z = otherVec.z;
};

vec3.prototype.add = function(otherVec)
{
        this.x += otherVec.x;
        this.y += otherVec.y;
        this.z += otherVec.z;
};

vec3.prototype.sub = function(otherVec)
{
        this.x -= otherVec.x;
        this.y -= otherVec.y;
        this.z -= otherVec.z;
};

vec3.prototype.mul = function(otherVec)
{
        this.x *= otherVec.x;
        this.y *= otherVec.y;
        this.z *= otherVec.z;
};

vec3.prototype.multScalar = function(scalarNumber)
{
        this.x *= scalarNumber;
        this.y *= scalarNumber;
        this.z *= scalarNumber;
};

vec3.prototype.squaredLength = function()
{
        return (this.x * this.x) + (this.y + this.y) + (this.z * this.z);
};

vec3.prototype.magnitude = function()
{
        return Math.hypot(this.x, this.y, this.z);
};

vec3.prototype.normalize = function()
{
        inverseMag = 1.0 / this.magnitude();
        this.x *= inverseMag;
        this.y *= inverseMag;
        this.z *= inverseMag;
};

vec3.prototype.distanceTo = function(otherVec)
{
        tempVec.copy(otherVec);
        tempVec.sub(this);
        return tempVec.magnitude();
};

vec3.prototype.dot = function(otherVec)
{
        return (this.x * otherVec.x) + (this.y * otherVec.y) + (this.z * otherVec.z);
};

vec3.prototype.crossVectors = function(vecA, vecB)
{
        this.x = (vecA.y * vecB.z) - (vecA.z * vecB.y);
        this.y = (vecA.z * vecB.x) - (vecA.x * vecB.z);
        this.z = (vecA.x * vecB.y) - (vecA.y * vecB.x);
};

vec3.prototype.randomCosWeightedDirectionInHemisphere = function(nl)
{
	let up = Math.sqrt(Math.random()); // cos-weighted distribution in hemisphere
    	let over = Math.sqrt(Math.max(0, 1 - (up * up)));
	let around = Math.random() * (Math.PI * 2);

	// from "Building an Orthonormal Basis, Revisited" http://jcgt.org/published/0006/01/01/
	let signf = nl.z >= 0 ? 1 : -1;
	let a = -1 / (signf + nl.z);
	let b = nl.x * nl.y * a;
	let T = new vec3( 1 + signf * nl.x * nl.x * a, signf * b, -signf * nl.x );
	let B = new vec3( b, signf + nl.y * nl.y * a, -nl.y );

        T.multScalar(Math.cos(around) * over);
        B.multScalar(Math.sin(around) * over);
        this.copy(nl);
        this.multScalar(up);
        this.add(T);
        this.add(B);
};

vec3.prototype.reflect = function(surfaceNormal)
{
        surfaceNormal.normalize();
        sN.copy(surfaceNormal);
        sN.multScalar(2 * this.dot(surfaceNormal));
        this.sub(sN);
};

let IdotN = 0.0;
let k = 0.0;
vec3.prototype.refract = function(surfaceNormal, eta)
{
        surfaceNormal.normalize();
        sN.copy(surfaceNormal);
        IdotN = this.dot(surfaceNormal);
        k = 1.0 - (eta * eta) * (1.0 - (IdotN * IdotN));
        sN.multScalar(eta * IdotN + Math.sqrt(k));
        this.multScalar(eta);
        this.sub(sN);
};

let temp = 0.0;
let cosi = 0.0;
let sint = 0.0;
let cost = 0.0;
let Rs = 0.0;
let Rp = 0.0;
function calcFresnelReflectance(rayDirection, n, etai, etat)
{
	temp = etai;
	cosi = rayDirection.dot(n);
	if (cosi > 0.0)
	{
		etai = etat;
		etat = temp;
	}
	
	ratioIoR = etai / etat;
	sint = ratioIoR * Math.sqrt(1.0 - (cosi * cosi));
	if (sint > 1.0) 
		return 1.0; // total internal reflection

	cost = Math.sqrt(1.0 - (sint * sint));
	cosi = Math.abs(cosi);
	Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
	Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));

	return ((Rs * Rs) + (Rp * Rp)) * 0.5;
}

function tentFilter(x)
{
	return (x < 0.5) ? Math.sqrt(2.0 * x) - 1.0 : 1.0 - Math.sqrt(2.0 - (2.0 * x));
}




function solveQuadratic(A, B, C)
{
	discrim = (B * B) - 4.0 * (A * C);
    
	if (discrim < 0.0)
        	return false;
    
	rootDiscrim = Math.sqrt(discrim);

	Q = (B > 0.0) ? -0.5 * (B + rootDiscrim) : -0.5 * (B - rootDiscrim); 
	
        t0 = C / Q;
        t1 = Q / A;
	
	return true;
}

function SphereIntersect( radius, position, rayOrigin, rayDirection )
{
        t0 = 0;
        t1 = 0;

        L.copy(rayOrigin);
	L.sub(position);
	a = rayDirection.dot(rayDirection);
	b = 2.0 * rayDirection.dot(L);
	c = L.dot(L) - (radius * radius);
        
        if ( !solveQuadratic(a, b, c) )
                return Infinity;

	if (t0 > 0.0)
		return t0;

	if (t1 > 0.0)
		return t1;
		
	return Infinity;
}

let denom = 0.0;
let pOrO = new vec3();
let scaledNormal = new vec3();
let result = 0.0;
function PlaneIntersect( planeNormal, d, rayOrigin, rayDirection )
{
	planeNormal.normalize();
	denom = planeNormal.dot(rayDirection);

	// uncomment the following if single-sided plane is desired
	//if (denom >= 0.0) return Infinity;
        
        scaledNormal.copy(planeNormal);
        scaledNormal.multScalar(d);
        pOrO.copy(scaledNormal);
        pOrO.sub(rayOrigin); 
        result = pOrO.dot(planeNormal) / denom;
	return (result > 0.0) ? result : Infinity;
}


let canvas = document.querySelector('canvas');
canvas.width = 640;//document.body.clientWidth;
canvas.height = 480;//document.body.clientHeight;
let canvasWidth = canvas.width;
let canvasHeight = canvas.height;
let invWidth = 1.0 / canvas.width;
let invHeight = 1.0 / canvas.height;
let aspectRatio = canvas.width / canvas.height;
let FOV = 60.0;
let thetaFOV = FOV * 0.5 * (Math.PI / 180.0);
let vLen = Math.tan(thetaFOV); // height scale
let uLen = vLen * aspectRatio; // width scale

let ctx = canvas.getContext('2d');
let imageData = ctx.getImageData(0, 0, canvas.width, canvas.height);

let pixelMemory = [];
for (let i = 0; i < imageData.data.length; i += 4) 
{
        pixelMemory[i] = new vec3();
}

let u = 0;
let v = 0;
let pixelOffsetX = 0;
let pixelOffsetY = 0;
let pixelColor = new vec3();
let colorPlusOne = new vec3();
let inverseColor = new vec3();
let cameraOrigin = new vec3(-2, 3.5, 6);//0,3,6
let cameraTarget = new vec3(0, 1, 0);
let cameraRightVec = new vec3();
let cameraUpVec = new vec3();
let cameraForwardVec = new vec3();
let worldUpVector = new vec3(0, 1, 0);
let tempUpVec = new vec3();
let tempRightVec = new vec3();
let rayOrigin = new vec3();
let rayDirection = new vec3();
let secondaryRayOrigin = new vec3();
let secondaryRayDirection = new vec3();
let accumulatedColor = new vec3();
let colorMask = new vec3();
let secondaryColorMask = new vec3();
let bounceIsSpecular = true;
let sampleLight = false;
let diffuseCount = 0;
let firstTypeWasDialectric = false;
let firstTypeWasDiffuse = false;
let isReflectionTime = false;
let isShadowTime = false;
let t = Infinity;
let d = Infinity;
let metalSphereRad = 2;
let metalSpherePos = new vec3(-2, metalSphereRad, -3);
let diffuseSphereRad = 1;
let diffuseSpherePos = new vec3(-2, diffuseSphereRad, 0.5);
let glassSphereRad = 1.5;
let glassSpherePos = new vec3(2, glassSphereRad, 1);
let coatSphereRad = 2;
let coatSpherePos = new vec3(3, coatSphereRad, -4);
let planeNormal = new vec3(0, 1, 0);
let planeD = 0;
let sunDirection = new vec3(1, 1, 0.5);
sunDirection.normalize();
let sunColor = new vec3(1.0, 0.95, 0.9);
sunColor.multScalar(5);
let skyColor = new vec3(0.3, 0.6, 1.0);
skyColor.multScalar(2);
let checkColor1 = new vec3(1, 1, 1);
let checkColor2 = new vec3(0, 0, 0);
let checkScale = 1.0;
let intersectionPoint = new vec3();
let hitPoint = new vec3();
let normal = new vec3();
let n = new vec3();
let nl = new vec3();
let nc = 1.0; // IOR of Air
let nt = 1.5; // IOR of common Glass
let Re = 0.0;
let Tr = 0.0;
let P = 0.0;
let tempNormal = new vec3();
let tempDir = new vec3();
let hitRecord = {};
hitRecord.t = Infinity;
hitRecord.type = -Infinity;
hitRecord.color = new vec3();
hitRecord.normal = new vec3();
const CHECKER = 0;
const DIFFUSE = 1;
const METAL = 2;
const DIALECTRIC = 3;
const CLEARCOAT = 4;
const MAX_BOUNCES = 6;

// each sample takes about a second, so use this number with caution!
const MAX_SAMPLE_COUNT = 50;// this number * 1 second = total rendering time to finish
let sampleCount = 0;
let infoElement = document.getElementById("info");


function sceneIntersect(rayOrigin, rayDirection)
{
        hitRecord.t = Infinity;
        d = Infinity;

        d = SphereIntersect(metalSphereRad, metalSpherePos, rayOrigin, rayDirection);
        if (d < hitRecord.t)
        {
                hitRecord.t = d;
                hitRecord.color.set(1.0, 0.8, 0.2);
                tempDir.copy(rayDirection);
                tempDir.multScalar(hitRecord.t);
                hitPoint.copy(rayOrigin);
                hitPoint.add(tempDir);
                hitRecord.normal.copy(hitPoint);
                hitRecord.normal.sub(metalSpherePos);
                hitRecord.normal.normalize();
                hitRecord.type = METAL;
        }

        d = SphereIntersect(diffuseSphereRad, diffuseSpherePos, rayOrigin, rayDirection);
        if (d < hitRecord.t)
        {
                hitRecord.t = d;
                hitRecord.color.set(1, 0, 0);
                tempDir.copy(rayDirection);
                tempDir.multScalar(hitRecord.t);
                hitPoint.copy(rayOrigin);
                hitPoint.add(tempDir);
                hitRecord.normal.copy(hitPoint);
                hitRecord.normal.sub(diffuseSpherePos);
                hitRecord.normal.normalize();
                hitRecord.type = DIFFUSE;
        }

        d = SphereIntersect(glassSphereRad, glassSpherePos, rayOrigin, rayDirection);
        if (d < hitRecord.t)
        {
                hitRecord.t = d;
                hitRecord.color.set(0.6, 1.0, 0.9);
                tempDir.copy(rayDirection);
                tempDir.multScalar(hitRecord.t);
                hitPoint.copy(rayOrigin);
                hitPoint.add(tempDir);
                hitRecord.normal.copy(hitPoint);
                hitRecord.normal.sub(glassSpherePos);
                hitRecord.normal.normalize();
                hitRecord.type = DIALECTRIC;
        }

        d = SphereIntersect(coatSphereRad, coatSpherePos, rayOrigin, rayDirection);
        if (d < hitRecord.t)
        {
                hitRecord.t = d;
                hitRecord.color.set(0.0, 0.0, 1.0);
                tempDir.copy(rayDirection);
                tempDir.multScalar(hitRecord.t);
                hitPoint.copy(rayOrigin);
                hitPoint.add(tempDir);
                hitRecord.normal.copy(hitPoint);
                hitRecord.normal.sub(coatSpherePos);
                hitRecord.normal.normalize();
                hitRecord.type = CLEARCOAT;
        }

        d = PlaneIntersect(planeNormal, planeD, rayOrigin, rayDirection);
        if (d < hitRecord.t)
        {
                hitRecord.t = d;
                hitRecord.color.set(1, 1, 1);
                tempDir.copy(rayDirection);
                tempDir.multScalar(hitRecord.t);
                //hitPoint.copy(rayOrigin);
                //hitPoint.add(tempDir);
                hitRecord.normal.copy(planeNormal);
                hitRecord.normal.normalize();
                hitRecord.type = CHECKER;
        }

        return hitRecord;
}

function rayTrace(rayOrigin, rayDirection)
{
        accumulatedColor.set(0, 0, 0);
        colorMask.set(1, 1, 1);
        diffuseCount = 0;
        bounceIsSpecular = true;
        sampleLight = false;
        firstTypeWasDialectric = false;
        firstTypeWasDiffuse = false;
        isReflectionTime = false;
        isShadowTime = false;

        for (let bounces = 0; bounces < MAX_BOUNCES; bounces++)
        {
                hitRecord = sceneIntersect(rayOrigin, rayDirection);

                if (hitRecord.t == Infinity)
                {
                        if (bounces == 0)
                        {
                                accumulatedColor.copy(colorMask);
                                accumulatedColor.mul(skyColor);
                                break;
                        }

                        if (firstTypeWasDialectric)
			{
				if (!isReflectionTime) 
				{
                                        if (bounceIsSpecular)
                                        {
                                                accumulatedColor.copy(colorMask);
                                                accumulatedColor.mul(skyColor);
                                        }

                                        if (sampleLight)
                                        {
                                                accumulatedColor.copy(colorMask);
                                                accumulatedColor.mul(sunColor);
                                        }
					
					// start back at the refractive surface, but this time follow reflective branch
                                        rayOrigin.copy(secondaryRayOrigin);
                                        rayDirection.copy(secondaryRayDirection);
					colorMask.copy(secondaryColorMask);
					// set/reset variables
					isReflectionTime = true;
					bounceIsSpecular = true;
					sampleLight = false;
					// continue with the reflection ray
					continue;
				} // end if (!isReflectionTime) 

                                // add reflective result to the refractive result (if any)
                                if (bounceIsSpecular)
                                {
                                        colorMask.mul(skyColor);
                                        accumulatedColor.add(colorMask);
                                        break;
                                }

                                if (sampleLight)
                                {
                                        colorMask.mul(sunColor);
                                        accumulatedColor.add(colorMask);
                                        break;
                                }
                                
			} // end if (firstTypeWasDialectric)


                        if (firstTypeWasDiffuse)
                        {
                                if (!isShadowTime)
                                {
                                        accumulatedColor.copy(colorMask);
                                        accumulatedColor.mul(skyColor);
                                        accumulatedColor.multScalar(0.5);

                                        // start back at the diffuse surface, but this time follow shadow ray branch
                                        rayOrigin.copy(secondaryRayOrigin);
                                        rayDirection.copy(secondaryRayDirection);
                                        colorMask.copy(secondaryColorMask);
                                        // set/reset variables
                                        isShadowTime = true;
                                        bounceIsSpecular = false;
                                        sampleLight = true;
                                        // continue with the shadow ray
                                        continue;
                                }
                                
                                colorMask.mul(sunColor);
                                colorMask.multScalar(0.5);
                                accumulatedColor.add(colorMask);
                                break;
                        } // end if (firstTypeWasDiffuse)
                        

                        if (bounceIsSpecular)
                        {
                                accumulatedColor.copy(colorMask);
                                accumulatedColor.mul(skyColor);
                                break;
                        }

                } // end if (hitRecord.t == Infinity)


                // if we reached this point and sampleLight is still true, that means the shadow ray ran into other scene
                // geometry before it could reach the light source, so continue with secondary rays, or if already doing that, exit 
                if (sampleLight)
                {
                        if (firstTypeWasDialectric && !isReflectionTime) 
			{
				// start back at the refractive surface, but this time follow reflective branch
				rayOrigin.copy(secondaryRayOrigin);
                                rayDirection.copy(secondaryRayDirection);
                                colorMask.copy(secondaryColorMask);
                                // set/reset variables
                                isReflectionTime = true;
                                bounceIsSpecular = true;
                                sampleLight = false;
                                // continue with the reflection ray
                                continue;
                        }
                        
                        if (firstTypeWasDiffuse && !isShadowTime) 
			{
				// start back at the diffuse surface, but this time follow shadow ray branch
				rayOrigin.copy(secondaryRayOrigin);
                                rayDirection.copy(secondaryRayDirection);
                                colorMask.copy(secondaryColorMask);
                                // set/reset variables
                                isShadowTime = true;
                                bounceIsSpecular = false;
                                sampleLight = true;
                                // continue with the shadow ray
                                continue;
			}

			// nothing left to calculate, so exit	
                        break; // this exit creates a shadow
                }

                // useful data 
                n.copy(hitRecord.normal);
                n.normalize();
                nl.copy(n);
                if (rayDirection.dot(n) >= 0.0)
                        nl.multScalar(-1);
                nl.normalize();

                // calculate intersection point
                tempDir.copy(rayDirection);
                tempDir.multScalar(hitRecord.t);
                rayOrigin.add(tempDir); // ray origin is now located at the intersection point
                tempNormal.copy(nl);
                tempNormal.multScalar(0.01);


                if (hitRecord.type == CHECKER)
                {
                        diffuseCount++;

                        // create checker pattern
                        //if ( Math.abs( Math.floor(rayOrigin.x * checkScale) ) % 2 + Math.abs( Math.floor(rayOrigin.z * checkScale) ) % 2 == 1 )
                        if ( Math.sin((rayOrigin.x * checkScale)) > -0.98 && Math.sin((rayOrigin.z * checkScale)) > -0.98 )    
                                hitRecord.color.copy(checkColor1);
                        else hitRecord.color.copy(checkColor2);

                        colorMask.mul(hitRecord.color);

                        bounceIsSpecular = false;

                        if (diffuseCount == 1 && !firstTypeWasDiffuse && !firstTypeWasDialectric)
			{	
				// save intersection data for future shadowray trace
                                firstTypeWasDiffuse = true;
                                secondaryRayOrigin.copy(rayOrigin);
                                secondaryRayOrigin.add(tempNormal);
                                secondaryRayDirection.copy(sunDirection); // create shadow ray
                                secondaryRayDirection.normalize();
                                secondaryColorMask.copy(colorMask);
                                secondaryColorMask.multScalar(Math.max(0.0, secondaryRayDirection.dot(nl)));

                                // choose random Diffuse sample vector
                                rayOrigin.add(tempNormal);
                                rayDirection.randomCosWeightedDirectionInHemisphere(nl);
                                rayDirection.normalize();
                                continue;
                        }

                        rayOrigin.add(tempNormal);
                        rayDirection.copy(sunDirection); // create shadow ray
                        rayDirection.normalize();
                        colorMask.multScalar(Math.max(0.0, rayDirection.dot(nl)));
                        sampleLight = true;

                        continue;
                }

                if (hitRecord.type == DIFFUSE)
                { 
                        diffuseCount++;

                        colorMask.mul(hitRecord.color);

                        bounceIsSpecular = false;

                        if (diffuseCount == 1 && !firstTypeWasDiffuse && !firstTypeWasDialectric)
			{	
				// save intersection data for future shadowray trace
                                firstTypeWasDiffuse = true;
                                secondaryRayOrigin.copy(rayOrigin);
                                secondaryRayOrigin.add(tempNormal);
                                secondaryRayDirection.copy(sunDirection); // create shadow ray
                                secondaryRayDirection.normalize();
                                secondaryColorMask.copy(colorMask);
                                secondaryColorMask.multScalar(Math.max(0.0, secondaryRayDirection.dot(nl)));

                                // choose random Diffuse sample vector
                                rayOrigin.add(tempNormal);
                                rayDirection.randomCosWeightedDirectionInHemisphere(nl);
                                rayDirection.normalize();
                                continue;
                        }

                        rayOrigin.add(tempNormal);
                        rayDirection.copy(sunDirection); // create shadow ray
                        rayDirection.normalize();
                        colorMask.multScalar(Math.max(0.0, rayDirection.dot(nl)));
                        sampleLight = true;

                        continue;
                }

                if (hitRecord.type == METAL)
                { 
                        colorMask.mul(hitRecord.color);

                        rayOrigin.add(tempNormal); // nudge ray out from surface a tiny bit along surface normal
                        // the above is needed to prevent self-intersection with previously intersected surface

                        rayDirection.reflect(nl); // create reflection ray
                        rayDirection.normalize();
                        //bounceIsSpecular = true;// turn on mirror reflecting caustics

                        continue; 
                }

                if (hitRecord.type == DIALECTRIC)
                { 
                        nc = 1.0; // IOR of Air
			nt = 1.5; // IOR of common Glass
                        Re = calcFresnelReflectance(rayDirection, n, nc, nt);
                        Re = Math.max(0.0, Math.min(1.0, Re));
                        Tr = 1.0 - Re;
                        //P = (Re + 0.5) * 0.5;

                        if (diffuseCount == 0 && !firstTypeWasDialectric)
                        {
                                firstTypeWasDialectric = true;
                                // create secondary REFLECTION ray for a future reflection trace
                                secondaryRayOrigin.copy(rayOrigin);
                                secondaryRayOrigin.add(tempNormal);
                                secondaryRayDirection.copy(rayDirection);
                                secondaryRayDirection.reflect(nl);
                                secondaryRayDirection.normalize();
                                secondaryColorMask.copy(colorMask);
                                secondaryColorMask.multScalar(Re);

                                colorMask.multScalar(Tr);
                        }
                        
                        // REFRACT (Transmit)
                        colorMask.mul(hitRecord.color);
                        
                        rayOrigin.sub(tempNormal);
                        rayDirection.refract(nl, ratioIoR);
                        rayDirection.normalize();

                        bounceIsSpecular = true;// turn on refracting caustics

                        continue; 
                } // end if (hitRecord.type == DIALECTRIC)

                if (hitRecord.type == CLEARCOAT)
                { 
                        nc = 1.0; // IOR of Air
			nt = 1.4; // IOR of clearCoat
                        Re = calcFresnelReflectance(rayDirection, n, nc, nt);
                        Re = Math.max(0.0, Math.min(1.0, Re));
                        Tr = 1.0 - Re;

                        // clearCoat counts as dialectric
                        if (diffuseCount == 0 && !firstTypeWasDialectric)
                        {
                                firstTypeWasDialectric = true;
                                // create secondary REFLECTION ray for a future reflection trace
                                secondaryRayOrigin.copy(rayOrigin);
                                secondaryRayOrigin.add(tempNormal);
                                secondaryRayDirection.copy(rayDirection);
                                secondaryRayDirection.reflect(nl);
                                secondaryRayDirection.normalize();
                                secondaryColorMask.copy(colorMask);
                                secondaryColorMask.multScalar(Re);

                                colorMask.multScalar(Tr);
                        }
                        
                        
                        diffuseCount++;

                        colorMask.mul(hitRecord.color);

                        //colorMask.multScalar(2);

                        bounceIsSpecular = false;
                        
                        if (bounces == 0 && Math.random() < 0.5)
                        {
                                // choose random Diffuse sample vector
                                rayOrigin.add(tempNormal);
                                rayDirection.randomCosWeightedDirectionInHemisphere(nl);
                                rayDirection.normalize();
                                continue;
                        }
                        else
                        {
                                // send shadow ray
                                rayOrigin.add(tempNormal);
                                rayDirection.copy(sunDirection);
                                rayDirection.normalize();
                                colorMask.multScalar(Math.max(0.0, rayDirection.dot(nl)));
                                sampleLight = true;
                                continue;
                        }
                        

                } // end if (hitRecord.type == DIALECTRIC)

        } // end for (let bounces = 0; bounces < MAX_BOUNCES; bounces++)

        return accumulatedColor;
} // end function rayTrace(rayOrigin, rayDirection)


function getPixelColor()
{
        // construct 'lookAt' camera frame
        cameraForwardVec.copy(cameraOrigin);
        cameraForwardVec.sub(cameraTarget);
        cameraForwardVec.normalize();
        cameraRightVec.crossVectors(worldUpVector, cameraForwardVec);
        cameraRightVec.normalize();
        cameraUpVec.crossVectors(cameraForwardVec, cameraRightVec);
        cameraUpVec.normalize();
        // the camera's forward vec needed to be constructed from the camera target pointing towards the camera location
        // for generating the mathematically correct orthonormal basis of the camera.  But the actual ray direction needs to
        // point in the opposite direction, in other words from the camera location pointing towards the camera target
        cameraForwardVec.multScalar(-1); // flip it to opposite direction
        cameraForwardVec.normalize();

        imageData = ctx.getImageData(0, 0, canvas.width, canvas.height);
        
        // loop over every pixel on the canvas all the way from the top left to the bottom right
        for (let i = 0; i < imageData.data.length; i += 4) 
        {      
 
                u = (((i * 0.25) % canvasWidth) * invWidth) * 2.0 - 1.0;
                // left side of image is -1, middle is 0, and right is +1 with smooth transition in between
                v = (Math.floor((i * 0.25) * invWidth) * invHeight) * 2.0 - 1.0;
                v *= -1.0;// flip Y(v) coordinates, bottom of image is now -1, middle is 0, and top is +1 with smooth transition

                pixelOffsetX = tentFilter(Math.random());
                pixelOffsetY = tentFilter(Math.random());
                pixelOffsetX /= (canvasWidth);
                pixelOffsetY /= (canvasHeight);

                u += pixelOffsetX;
                v += pixelOffsetY;

                // construct ray for this particular pixel (u,v)
                rayOrigin.copy(cameraOrigin);

                tempRightVec.copy(cameraRightVec);
                tempRightVec.multScalar(u);
                tempRightVec.multScalar(uLen);
                tempUpVec.copy(cameraUpVec);
                tempUpVec.multScalar(v);
                tempUpVec.multScalar(vLen);
                
                rayDirection.copy(cameraForwardVec);
                rayDirection.add(tempRightVec);
                rayDirection.add(tempUpVec);
                rayDirection.normalize();

                // calculate this pixel's color through ray tracing
                pixelMemory[i].add(rayTrace(rayOrigin, rayDirection));
                pixelColor.copy(pixelMemory[i]);

                pixelColor.multScalar(1.0 / (sampleCount + 1));

                // apply Reinhard tonemapping
                colorPlusOne.set(1.0 + pixelColor.x, 1.0 + pixelColor.y, 1.0 + pixelColor.z);
                inverseColor.set(1.0 / colorPlusOne.x, 1.0 / colorPlusOne.y, 1.0 / colorPlusOne.z);
                pixelColor.mul(inverseColor);

                // clamp values to 0-1 range
                // pixelColor.x = Math.max(0.0, Math.min(1.0, pixelColor.x));
                // pixelColor.y = Math.max(0.0, Math.min(1.0, pixelColor.y));
                // pixelColor.z = Math.max(0.0, Math.min(1.0, pixelColor.z));
	
                // do gamma correction
                pixelColor.set(Math.pow(pixelColor.x, 0.4545), Math.pow(pixelColor.y, 0.4545), Math.pow(pixelColor.z, 0.4545)); 

                
                imageData.data[i + 0] = pixelColor.x * 255; // red
                imageData.data[i + 1] = pixelColor.y * 255; // green
                imageData.data[i + 2] = pixelColor.z * 255; // blue
                imageData.data[i + 3] = 255;                // alpha
        }

        ctx.putImageData(imageData, 0, 0);
        sampleCount++;
}

function animate() 
{
        getPixelColor();
        infoElement.innerHTML = "Samples: " + sampleCount + " / " + MAX_SAMPLE_COUNT;

        if (sampleCount < MAX_SAMPLE_COUNT)
                requestAnimationFrame(animate);      
}

// start up the progressive rendering!
animate();
