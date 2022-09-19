/* 3D Vector/Surface Math library... */

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

vec3.prototype.negate = function ()
{
	this.x *= -1;
	this.y *= -1;
	this.z *= -1;
};

vec3.prototype.mix = function(vecA, vecB, amount)
{
	amount = Math.max(0, amount); // clamp supplied amount to 0-1 range
	amount = Math.min(1, amount); // clamp supplied amount to 0-1 range
	this.x = (vecA.x * (1 - amount)) + (vecB.x * amount);
	this.y = (vecA.y * (1 - amount)) + (vecB.y * amount);
	this.z = (vecA.z * (1 - amount)) + (vecB.z * amount);
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

vec3.prototype.reflect = function(surfaceNormal)
{
	// R = I + (N * 2 * -IdotN)
	//surfaceNormal.normalize();
	sN.copy(surfaceNormal);
	sN.multScalar(2 * -this.dot(surfaceNormal));
	this.add(sN);
};

let IdotN = 0.0;
let k = 0.0;
vec3.prototype.refract = function(surfaceNormal, eta)
{
	// T = eta * I + (N * eta * -IdotN - sqrt(1 - (eta * eta * (1 - IdotN * IdotN))))
	//surfaceNormal.normalize();
	sN.copy(surfaceNormal);
	IdotN = -this.dot(surfaceNormal);
	k = (eta * eta) * (1 - IdotN * IdotN);
	sN.multScalar(eta * IdotN - Math.sqrt(1 - k));
	this.multScalar(eta);
	this.add(sN);
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

let halfDir = new vec3();
let specAngle = 0;
let specAmount = 0;
let specularContribution = new vec3();

function calcBlinnPhongReflection(rayDirection, normal, shininessExp)
{
	halfDir.copy(rayDirection);
	halfDir.negate();
	halfDir.add(sunDirection);
	halfDir.normalize();
	
	specAngle = Math.max(0.0, halfDir.dot(normal));
	
	return Math.pow(specAngle, shininessExp);
}


function tentFilter(x)
{
	return (x < 0.5) ? Math.sqrt(2.0 * x) - 1.0 : 1.0 - Math.sqrt(2.0 - (2.0 * x));
}


let invA, neg_halfB, u2, u1;
function solveQuadratic(A, B, C)
{
	invA = 1.0 / A;
	B *= invA;
	C *= invA;
	neg_halfB = -B * 0.5;
	u2 = neg_halfB * neg_halfB - C;
	u1 = u2 < 0.0 ? neg_halfB = 0.0 : Math.sqrt(u2);
	t0 = neg_halfB - u1;
	t1 = neg_halfB + u1;
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
	//planeNormal.normalize();
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

let numOfDownloads = 0;
let canvas = document.querySelector('canvas');
document.querySelector('a').addEventListener('click', function(event) {
	event.target.href = canvas.toDataURL();
	event.target.download = "render" + numOfDownloads + ".png";
	numOfDownloads++;	
});
canvas.width = window.innerWidth;
canvas.height = window.innerHeight;
let invWidth = 1.0 / canvas.width;
let invHeight = 1.0 / canvas.height;
let aspectRatio = canvas.width / canvas.height;
let FOV = 30.0;
let thetaFOV = FOV * 0.5 * (Math.PI / 180.0);
let vLen = Math.tan(thetaFOV); // height scale
let uLen = vLen * aspectRatio; // width scale

let ctx = canvas.getContext('2d');
let imageData = ctx.getImageData(0, 0, canvas.width, canvas.height);

let RayTracedColorHistory = [];

for (let j = 0; j < canvas.height; j++) 
{
	for (let i = 0; i < canvas.width; i++) 
	{
		RayTracedColorHistory[(j * canvas.width) + i] = new vec3();
	}
}


let u = 0;
let v = 0;
let pixelOffsetX = 0;
let pixelOffsetY = 0;
let pixelIndex = 0;
let pixelColor = new vec3();
let colorPlusOne = new vec3();
let inverseColor = new vec3();

//let cameraOrigin = new vec3(2, 2, 1.5);
//let cameraTarget = new vec3(1, 1.9, 1);
let cameraOrigin = new vec3(-1.2, 0.7, 3);
cameraOrigin.normalize();
cameraOrigin.multScalar(15);
let cameraTarget = new vec3(0, 1, 0);

let cameraRightVec = new vec3();
let cameraUpVec = new vec3();
let cameraForwardVec = new vec3();
let worldUpVector = new vec3(0, 1, 0);
let tempUpVec = new vec3();
let tempRightVec = new vec3();
let rayOrigin = new vec3();
let rayDirection = new vec3();
let ambientColor = new vec3();
let diffuseColor = new vec3();
let specularColor = new vec3();
let accumulatedColor = new vec3();
let colorMask = new vec3();
let bounceIsSpecular = true;
let sendShadowRay = false;
let t = Infinity;
let d = Infinity;
let metalSphereRad = 2;
let metalSpherePos = new vec3(-2, metalSphereRad, -3);
let diffuseSphereRad = 1;
let diffuseSpherePos = new vec3(-2, diffuseSphereRad, 2);
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
let checkColor1 = new vec3(0.5, 0.5, 0.0);
let checkColor2 = new vec3(0.5, 0.0, 0.0);
let checkScale = 0.3;
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
let RP = 0.0;
let TP = 0.0;
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
const TRANSPARENT = 3;
const CLEARCOAT = 4;
const MAX_BOUNCES = 6;

// each sample takes about a second, so use this number with caution!
const MAX_SAMPLE_COUNT = 100;// this number * 1 second = total rendering time to finish
let sampleCount = 0;
let infoElement = document.getElementById("info");

window.addEventListener('resize', onWindowResize, false);

function onWindowResize(event)
{
	// reset camera
	cameraOrigin.set(-1.2, 0.7, 3);
	cameraOrigin.normalize();
	cameraOrigin.multScalar(15);
	// recalculate image dimensions and viewing plane scale
	canvas.width = window.innerWidth;
	canvas.height = window.innerHeight;
	invWidth = 1.0 / canvas.width;
	invHeight = 1.0 / canvas.height;
	aspectRatio = canvas.width / canvas.height;
	vLen = Math.tan(thetaFOV); // height scale
	uLen = vLen * aspectRatio; // width scale

	imageData = ctx.getImageData(0, 0, canvas.width, canvas.height);

	RayTracedColorHistory = [];

	for (let j = 0; j < canvas.height; j++) 
	{
		for (let i = 0; i < canvas.width; i++) 
		{
			RayTracedColorHistory[(j * canvas.width) + i] = new vec3();
		}
	}

	if (sampleCount == MAX_SAMPLE_COUNT)
	{
		sampleCount = 0; // reset sampleCount to start a fresh progressive render
		animate(); // also need to restart requestAnimationFrame
	}
	else
		sampleCount = 0; // reset sampleCount to start a fresh progressive render
}

function sceneIntersect(rayOrigin, rayDirection)
{
	hitRecord.t = Infinity;
	d = Infinity;

	d = SphereIntersect(metalSphereRad, metalSpherePos, rayOrigin, rayDirection);
	if (d < hitRecord.t)
	{
		hitRecord.t = d;
		hitRecord.color.set(1.0, 0.7, 0.1);
		tempDir.copy(rayDirection);
		tempDir.multScalar(hitRecord.t);
		hitPoint.copy(rayOrigin);
		hitPoint.add(tempDir);
		hitRecord.normal.copy(hitPoint);
		hitRecord.normal.sub(metalSpherePos);
		hitRecord.type = METAL;
	}

	d = SphereIntersect(diffuseSphereRad, diffuseSpherePos, rayOrigin, rayDirection);
	if (d < hitRecord.t)
	{
		hitRecord.t = d;
		hitRecord.color.set(1.0, 0.0, 0.0);
		tempDir.copy(rayDirection);
		tempDir.multScalar(hitRecord.t);
		hitPoint.copy(rayOrigin);
		hitPoint.add(tempDir);
		hitRecord.normal.copy(hitPoint);
		hitRecord.normal.sub(diffuseSpherePos);
		hitRecord.type = DIFFUSE;
	}

	d = SphereIntersect(glassSphereRad, glassSpherePos, rayOrigin, rayDirection);
	if (d < hitRecord.t)
	{
		hitRecord.t = d;
		//hitRecord.color.set(0.6, 1.0, 0.9);
		hitRecord.color.set(0.2, 1.0, 1.0);
		tempDir.copy(rayDirection);
		tempDir.multScalar(hitRecord.t);
		hitPoint.copy(rayOrigin);
		hitPoint.add(tempDir);
		hitRecord.normal.copy(hitPoint);
		hitRecord.normal.sub(glassSpherePos);
		hitRecord.type = TRANSPARENT;
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
		hitRecord.type = CLEARCOAT;
	}

	d = PlaneIntersect(planeNormal, planeD, rayOrigin, rayDirection);
	if (d < hitRecord.t)
	{
		hitRecord.t = d;
		hitRecord.color.set(1.0, 1.0, 1.0);
		tempDir.copy(rayDirection);
		tempDir.multScalar(hitRecord.t);
		//hitPoint.copy(rayOrigin);
		//hitPoint.add(tempDir);
		hitRecord.normal.copy(planeNormal);
		hitRecord.type = CHECKER;
	}

	return hitRecord;
} // end function sceneIntersect(rayOrigin, rayDirection)


function rayTrace(rayOrigin, rayDirection)
{
	accumulatedColor.set(0, 0, 0);
	ambientColor.set(0, 0, 0);
	diffuseColor.set(0, 0, 0);
	specularColor.set(0, 0, 0);
	specularContribution.set(0, 0, 0);
	colorMask.set(1, 1, 1);
	bounceIsSpecular = true;
	sendShadowRay = false;

	for (let bounces = 0; bounces < MAX_BOUNCES; bounces++)
	{
		hitRecord = sceneIntersect(rayOrigin, rayDirection);

		if (hitRecord.t == Infinity)
		{
			if (bounceIsSpecular)
			{
				accumulatedColor.copy(colorMask);
				accumulatedColor.mul(skyColor);
				accumulatedColor.add(specularContribution);
			}
			if ( !bounceIsSpecular )
			{
				accumulatedColor.add(specularContribution);
			}
			if (bounces == 0)
			{
				accumulatedColor.copy(skyColor);
			}
			
			break;        
		}


		// if we reached this point and sendShadowRay is still true, this means that the shadow ray 
		//  intersected another scene object before it could reach the light source, so exit
		if (sendShadowRay)
		{	
			accumulatedColor.copy(ambientColor);
			break; // this exit leaves a shadow
		}

		// useful data 
		n.copy(hitRecord.normal);
		n.normalize(); // this normalization is important! Hit normals come in from sceneIntersect() as non-normalized
		nl.copy(n);
		if (rayDirection.dot(n) >= 0.0)
			nl.multScalar(-1);

		// calculate intersection point
		tempDir.copy(rayDirection);
		tempDir.multScalar(hitRecord.t);
		rayOrigin.add(tempDir); // ray origin is now located at the intersection point
		tempNormal.copy(nl);
		tempNormal.multScalar(0.01);


		if (hitRecord.type == DIFFUSE || hitRecord.type == CHECKER)
		{
			if (hitRecord.type == CHECKER)
			{
				// create checker pattern
				if ( Math.abs( Math.floor(rayOrigin.x * checkScale) ) % 2 + Math.abs( Math.floor(rayOrigin.z * checkScale) ) % 2 == 1 )
				////if (Math.sin((rayOrigin.x * checkScale)) > -0.99 && Math.sin((rayOrigin.z * checkScale)) > -0.99)
					hitRecord.color.copy(checkColor1);
				else hitRecord.color.copy(checkColor2);
			}
			
			colorMask.mul(hitRecord.color);
			
			// ambient
			ambientColor.copy(colorMask);
			ambientColor.mul(skyColor);
			ambientColor.multScalar(0.3);

			// diffuse
			diffuseColor.copy(colorMask);
			diffuseColor.mul(sunColor);

			accumulatedColor.mix(ambientColor, diffuseColor, Math.max(0, nl.dot(sunDirection)));

			// create shadow ray
			rayOrigin.add(tempNormal);
			rayDirection.copy(sunDirection);

			bounceIsSpecular = false;
			sendShadowRay = true;
			continue;
		} // end if (hitRecord.type == DIFFUSE || hitRecord.type == CHECKER)


		if (hitRecord.type == METAL)
		{ 
			colorMask.mul(hitRecord.color);

			// specular contribution
			specularColor.copy(sunColor);
			specularColor.mul(colorMask);

			// evaluate Blinn-Phong reflection model at this surface point
			specAmount = calcBlinnPhongReflection(rayDirection, nl, 1000);
			specularColor.multScalar(specAmount);
			specularContribution.add(specularColor);

			rayOrigin.add(tempNormal); // nudge ray out from surface a tiny bit along surface normal
			// the above is needed to prevent self-intersection with previously intersected surface

			rayDirection.reflect(nl); // create reflection ray
			rayDirection.normalize();
			continue; 
		} // end if (hitRecord.type == METAL)

		if (hitRecord.type == TRANSPARENT)
		{ 
			nc = 1.0; // IOR of Air
			nt = 1.5; // IOR of common Glass
			Re = calcFresnelReflectance(rayDirection, n, nc, nt);
			Tr = 1.0 - Re;
			P = 0.25 + (0.5 * Re);
			RP = Re / P;
			TP = Tr / (1.0 - P);

			if (Math.random() < P)
			{
				colorMask.multScalar(RP);

				// specular contribution
				specularColor.copy(sunColor);

				// evaluate Blinn-Phong reflection model at this surface point
				specAmount = calcBlinnPhongReflection(rayDirection, nl, 1000);
				specularColor.multScalar(specAmount);
				specularContribution.add(specularColor);

				rayOrigin.add(tempNormal);
				rayDirection.reflect(nl);
				rayDirection.normalize();
				continue;
			}
			
			// REFRACT (Transmit)
			colorMask.multScalar(TP);
			colorMask.mul(hitRecord.color);

			// specular contribution
			specularColor.copy(sunColor);
			specularColor.mul(colorMask);

			// evaluate Blinn-Phong reflection model at this surface point
			specAmount = calcBlinnPhongReflection(rayDirection, nl, 1000);
			specularColor.multScalar(specAmount);
			specularContribution.add(specularColor);

			rayOrigin.sub(tempNormal);
			rayDirection.refract(nl, ratioIoR);
			rayDirection.normalize();
			continue; 
		} // end if (hitRecord.type == TRANSPARENT)

		if (hitRecord.type == CLEARCOAT)
		{ 
			nc = 1.0; // IOR of Air
			nt = 1.4; // IOR of clearCoat
			Re = calcFresnelReflectance(rayDirection, nl, nc, nt);
			Tr = 1.0 - Re;
			P = 0.25 + (0.5 * Re);
			RP = Re / P;
			TP = Tr / (1.0 - P);

			if (Math.random() < P)
			{
				colorMask.multScalar(RP);

				// specular contribution
				specularColor.copy(sunColor);

				// evaluate Blinn-Phong reflection model at this surface point
				specAmount = calcBlinnPhongReflection(rayDirection, nl, 1000);
				specularColor.multScalar(specAmount);
				specularContribution.add(specularColor);

				rayOrigin.add(tempNormal);
				rayDirection.reflect(nl);
				rayDirection.normalize();
				continue;
			}

			colorMask.multScalar(TP);

			colorMask.mul(hitRecord.color);
			
			// ambient
			ambientColor.copy(colorMask);
			ambientColor.mul(skyColor);
			ambientColor.multScalar(0.3);

			// diffuse
			diffuseColor.copy(colorMask);
			diffuseColor.mul(sunColor);

			accumulatedColor.mix(ambientColor, diffuseColor, Math.max(0, nl.dot(sunDirection)));

			
			// create shadow ray
			rayOrigin.add(tempNormal);
			rayDirection.copy(sunDirection);

			bounceIsSpecular = false;
			sendShadowRay = true;
			continue;
		} // end if (hitRecord.type == CLEARCOAT)

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
	

	imageData = ctx.getImageData(0, 0, canvas.width, canvas.height);
	
	// loop over every pixel on the canvas all the way from the top left to the bottom right
	for (let j = 0; j < canvas.height; j++)
	{
		for (let i = 0; i < canvas.width; i++) 
		{      
			// calculate random pixel offset (anti-aliasing)
			pixelOffsetX = tentFilter(Math.random());
			pixelOffsetY = tentFilter(Math.random());

			u = (i + pixelOffsetX) / canvas.width;
			v = (j + pixelOffsetY) / canvas.height;

			u = (u * 2) - 1;
			v = (v * 2) - 1;
			v *= -1;// flip Y(v) coordinates
			// now, bottom of image is -1, middle is 0, and top is +1 with smooth transition
			
			// construct ray that starts at camera origin and goes right through this particular pixel (u,v) on the view plane
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

			// get the flat 1d index for this pixel from its initial 2d indices in screen coordinates (i,j)
			pixelIndex = (j * canvas.width) + i; // convert 2d indices to flat 1d index

			// calculate this pixel's color through ray tracing
			// add current raytraced pixel color result into accumulation buffer (color history) for this particular pixel
			RayTracedColorHistory[pixelIndex].add(rayTrace(rayOrigin, rayDirection));
		// make a copy of the current raytraced pixel color, so we can apply post processing to it below safely,
		// without polluting the original pixel color value above that is added to the accumulation buffer (color history) each frame
			pixelColor.copy(RayTracedColorHistory[pixelIndex]);
			// now we can safely play around with 'pixelColor' and display it to the screen
			pixelColor.multScalar(1.0 / (sampleCount + 1));

			// apply Reinhard tonemapping (brings unbounded linear color values into 0-1 range)
			colorPlusOne.set(1.0 + pixelColor.x, 1.0 + pixelColor.y, 1.0 + pixelColor.z);
			inverseColor.set(1.0 / colorPlusOne.x, 1.0 / colorPlusOne.y, 1.0 / colorPlusOne.z);
			pixelColor.mul(inverseColor);

			// do gamma correction
			pixelColor.set(Math.pow(pixelColor.x, 0.4545), Math.pow(pixelColor.y, 0.4545), Math.pow(pixelColor.z, 0.4545)); 

			/* 
			// test screen uv pattern
			u = u * 0.5 + 0.5;
			v = v * 0.5 + 0.5;
			imageData.data[pixelIndex + 0] = u * Math.abs(Math.sin(sampleCount * 0.1)) * 255;// pixelColor.x * 255; // red
			imageData.data[pixelIndex + 1] = u * Math.abs(Math.cos(sampleCount * 0.1)) * 255;// pixelColor.y * 255; // green
			imageData.data[pixelIndex + 2] = v * Math.abs(Math.sin((sampleCount + 5) * 0.05)) * 255// pixelColor.z * 255; // blue
			imageData.data[pixelIndex + 3] = 255; // alpha 
			*/

			// clamp values to 0-1 range
			pixelColor.x = Math.max(0.0, Math.min(1.0, pixelColor.x));
			pixelColor.y = Math.max(0.0, Math.min(1.0, pixelColor.y));
			pixelColor.z = Math.max(0.0, Math.min(1.0, pixelColor.z));

			/* 
			'imageData' manipulation method (faster) for setting each screen pixel's color, one at a time 
			*/
			imageData.data[pixelIndex * 4 + 0] = pixelColor.x * 255; // red
			imageData.data[pixelIndex * 4 + 1] = pixelColor.y * 255; // green
			imageData.data[pixelIndex * 4 + 2] = pixelColor.z * 255; // blue
			imageData.data[pixelIndex * 4 + 3] = 255;                // alpha
			
			/* 
			alt 'fillRect' method (2x slower than above) of painting tiny pixel-sized rectangles one at a time, 
			all the way from top-left corner to bottom-right corner
			note: in order to use this method, comment out the line "ctx.putImageData(imageData, 0, 0);" 10 lines below
			*/
			// pixelColor.multScalar(255);
			// ctx.fillStyle = "rgb(" + pixelColor.x + "," + pixelColor.y + "," + pixelColor.z + ")";
			// ctx.fillRect(i, j, 1, 1);

		} // end for (let i = 0; i < canvas.width; i++) 

	} // end for (let j = 0; j < canvas.height; j++)
	
	ctx.putImageData(imageData, 0, 0);
	sampleCount++;
} // end function getPixelColor()


let movementVec = new vec3(0,0,-1);
function animate() 
{
	getPixelColor();
	infoElement.innerHTML = "Samples: " + sampleCount + " / " + MAX_SAMPLE_COUNT;

	//cameraOrigin.add(movementVec);

	if (sampleCount < MAX_SAMPLE_COUNT)
		requestAnimationFrame(animate);      
}

// start up the progressive rendering!
animate();
